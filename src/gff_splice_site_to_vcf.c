#include <htslib/vcf.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>

#include "bcftools/gff.h"
#include "bcftools/regidx.h"

#include "logging/log.h"

#define SPLICE_SITE_MALLOC_START_COUNT 1000000
#define BASES "ACGT"

typedef enum {
    ACCEPTOR,
    DONOR
} SpliceSiteType;

typedef struct {
    const char *chr;
    uint32_t rid; // Chromosome/region encoding, for quick comparisons
    uint64_t pos;
    SpliceSiteType type;

    const char *gene_name;
    uint64_t gene_beg;
    uint64_t gene_end;
    uint32_t tid;
    int strand;
} SpliceSite;

typedef struct {
    uint32_t n, m;
    SpliceSite *a;
} SpliceSites;

typedef struct
{
    int64_t start, end;
} Range;

const char *chromosomes[24] = {
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY"
};

SpliceSite init_splice_site(const char *chr, const uint32_t rid, const uint64_t pos, const SpliceSiteType type, const char *gene_name, const uint64_t beg, const uint64_t end, const uint32_t tid, const int strand) {
    return (SpliceSite) { chr, rid, pos, type, gene_name, beg, end, tid, strand };
}

void splice_sites_init(const uint32_t m, SpliceSites *sites) {
    sites->n = 0;
    sites->m = m;
    sites->a = malloc(m * sizeof(SpliceSite));
}

void splice_sites_destroy(SpliceSites *sites) {
    free(sites->a);
    sites->a = NULL;
    sites->n = 0;
    sites->m = 0;
}

int find_next_site_with_type(const SpliceSites sites, const uint32_t tid, const uint32_t start, const uint32_t end, const int direction, const SpliceSiteType target_type, SpliceSite *out) {
    if (direction > 0) {
        for (uint i = start+1; i < end; i++) {
            SpliceSite site = sites.a[i];
            if (site.tid != tid) break;
            if (site.type != target_type) continue;

            *out = site;
            return 0;
        }
    } else {
        for (uint32_t i = start-1; i >= end; i--) {
            SpliceSite site = sites.a[i];
            if (site.tid != tid) break;
            if (site.type != target_type) continue;

            *out = site;
            return 0;
        }
    }

    return 1;
}

SpliceSites get_splice_sites_from_gff(const gff_t *gff) {
    SpliceSites sites;
    splice_sites_init(SPLICE_SITE_MALLOC_START_COUNT, &sites);

    regidx_t *transcripts = gff_get((gff_t *) gff, idx_tscript); // Removes const qualification for this call as it's not typed const, but it is a simple getter without consequences
    regitr_t *itr = regitr_init(transcripts);
    while (regitr_loop(itr)) {
        const gf_tscript_t *tr = regitr_payload(itr, gf_tscript_t *);

        // Only want protein coding and stranded transcripts
        if ((tr->type != GF_PROTEIN_CODING) | (tr->gene->type != GF_PROTEIN_CODING)) continue;
        if (tr->strand == STRAND_UNK) log_error("Transcript has an UNKNOWN strand. Skipping...", itr->seq);

        const char *chr = itr->seq;
        const uint32_t rid = tr->gene->iseq;

        // Only want things located on chromosomes
        if (strncmp(chr, "chr", 3) != 0) continue; // WARN: This is based on RUMC GFF files, they may not always be prefixed with chr
        if (sites.n + 2 >= sites.m) sites.a = realloc(sites.a, (sites.m *= 2) * sizeof(SpliceSite));

        // Extract stuff from pointers, no practical benefit, just looks a bit cleaner
        const char *gene_name = tr->gene->name;
        const uint64_t tr_beg = tr->gene->beg, tr_end = tr->gene->end;
        const uint32_t tid = tr->id;
        const int strand = tr->strand;

        for (int i = 0; i < tr->ncds; i++) {
            const gf_cds_t *cds = tr->cds[i];
            const uint64_t cds_beg = cds->beg;
            const uint64_t cds_end = cds->beg + cds->len - 1; // Offset by -1 to get closed end coordinate

            if (tr->strand == STRAND_FWD) {
                if (i != 0) sites.a[sites.n++] = init_splice_site(chr, rid, cds_beg, ACCEPTOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the first exon
                if (i != tr->ncds-1) sites.a[sites.n++] = init_splice_site(chr, rid, cds_end, DONOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the last exon
            } else if (tr->strand == STRAND_REV) {
                if (i != 0) sites.a[sites.n++] = init_splice_site(chr, rid, cds_beg, DONOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the last exon (approaches from 3' -> 5')
                if (i != tr->ncds-1) sites.a[sites.n++] = init_splice_site(chr, rid, cds_end, ACCEPTOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the first exon (approaches from 3' -> 5')
            }
        }

    }

    return sites;
}

gff_t *read_gff(const char *gff_path) {
    log_debug("Loading GFF from: %s", gff_path);
    gff_t *gff = gff_init(gff_path);
    gff_parse(gff);
    return gff;
}


int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Run as: ./gff_splice_sites_to_vcf <human_fa> <gff> <out.vcf>");
        exit(1);
    }

    const char *fasta_in_path = argv[1];
    log_debug("Loading FASTA from: %s", fasta_in_path);
    faidx_t *fai = fai_load(fasta_in_path);

    const char *gff_in_path = argv[2];
    gff_t *gff = read_gff(gff_in_path);
    log_debug("Loading GFF from: %s", gff_in_path);
    SpliceSites sites = get_splice_sites_from_gff(gff);
    log_debug("Found %i splice sites", sites.n);

    const char *vcf_out_path = argv[3];
    log_debug("Opening VCF output file at: %s", vcf_out_path);
    htsFile *fp = hts_open(vcf_out_path, "w");
    if (!fp) {
        log_error("Failed to open VCF output file");
        return 1;
    }

    bcf_hdr_t *hdr = bcf_hdr_init("w");
    if (!hdr) {
        log_error("Failed to initialize VCF output file header");
        hts_close(fp);
        return 1;
    }

    bcf_hdr_append(hdr, "##fileformat=VCFv4.2");
    for (int i = 0; i < 24; i++) {
        char contig_header[256];
        snprintf(contig_header, sizeof(contig_header), "##contig=<ID=%s>", chromosomes[i]);
        bcf_hdr_append(hdr, contig_header);
    }

    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    bcf_hdr_append(hdr, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

    if (bcf_hdr_add_sample(hdr, "DNAXX_GENIGMA") < 0) {
        log_error("Failed to add sample to VCF header");
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        return 1;
    }

    if (bcf_hdr_write(fp, hdr) < 0) {
        log_error("Failed to write VCF header");
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        return 1;
    }

    for (int i = 0; i < sites.n; i++) {
        SpliceSite site = sites.a[i];

        Range splice_site_range;
        if (site.type == DONOR) {
            if (site.strand == STRAND_FWD) {
                splice_site_range = (Range) { (int64_t) site.pos + 1, (int64_t) site.pos + 3 };
            } else { // site.strand == STRAND_REV
                splice_site_range = (Range) { (int64_t) site.pos - 2, (int64_t) site.pos };
            }
        } else { // site.type == ACCEPTOR
            if (site.strand == STRAND_FWD) {
                splice_site_range = (Range) { (int64_t) site.pos - 2, (int64_t) site.pos };
            } else { // site.strand == STRAND_REV
                splice_site_range = (Range) { (int64_t) site.pos + 1, (int64_t) site.pos + 3 };
            }
        }

        int len = 0;
        const char *ref = faidx_fetch_seq(fai, site.chr, splice_site_range.start, splice_site_range.end, &len);

        for (int p = 0; p < 2; p++) {
            const char ref_base = ref[p];
            const uint64_t position = splice_site_range.start + p;

            for (int b = 0; b < 4; b++) {
                if (ref_base == BASES[b]) continue;

                bcf1_t *v = bcf_init();
                const char alt_base = BASES[b];

                v->rid = bcf_hdr_id2int(hdr, BCF_DT_CTG, site.chr);
                v->pos = position;

                v->n_allele = 2;
                v->rlen = 1;

                bcf_unpack(v, BCF_UN_ALL);

                v->n_sample = 1;
                v->d.als = malloc(4);
                v->d.m_als = 4;
                v->d.als[0] = ref_base;
                v->d.als[1] = '\0';
                v->d.als[2] = alt_base;
                v->d.als[3] = '\0';

                v->d.allele = malloc(2 * sizeof(char *));
                v->d.m_allele = 2;
                v->d.allele[0] = &(v->d.als[0]);
                v->d.allele[1] = &(v->d.als[2]);

                int32_t gt_arr[2] = {bcf_gt_unphased(0), bcf_gt_phased(1)}; // 0|1 genotype

                bcf_update_genotypes(hdr, v, gt_arr, 2);

                bcf_write(fp, hdr, v);

                bcf_destroy(v);
            }
        }
    }

    hts_close(fp);
}

