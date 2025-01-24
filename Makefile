
PHASED_HGDP_GVCF_PATHS='gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr"*".filtered.SNV_INDEL.phased.shapeit5.bcf"*"'
bcftools?=~/bcftools/bin/bcftools
samtools?=~/samtools/bin/samtools

# requires `uv`
.PHONY: sync
sync:
	uv sync

.PHONY: download-data
download-data:
	@echo "Downloading data..."
	@mkdir -p data/bcf
	gsutil -m cp ./data/bcf

.PHONY: bcf-to-vcf
bcf-to-vcf:
	@echo "Converting BCF to VCF..."
	@mkdir -p data/vcf
	for bcf in data/bcf/*.bcf; do \
  		echo Converting $$bcf to VCF...; \
		$(bcftools) convert -O z $$bcf -o data/vcf/$$(basename $$bcf .bcf).vcf.bgz; \
	done

.PHONY: install-autoscaling-policies
install-autoscaling-policies:
	curl https://raw.githubusercontent.com/hail-is/dataproc-autoscaling-configs/master/add-all.sh | bash

# assumes hailctl installation
# assumes autoscaling policies added
.PHONY: start-cluster
start-cluster:
	@sh ./clusters/create.sh

gnomad_af_base=gnomad_af_out.ht
gnomad_pop_base=gnomad_pop_out.ht
gnomad_af_bucket_path=$(GCP_DIVREF_BUCKET)/$(gnomad_af_base)
gnomad_pop_bucket_path=$(GCP_DIVREF_BUCKET)/$(gnomad_pop_base)

.PHONY: extract-gnomad-afs
extract-gnomad-afs:
	@echo "Extracting gnomAD allele frequencies..."
	uv run hailctl dataproc submit haplo1 scripts/extract_gnomad_afs.py -- \
		$(gnomad_af_bucket_path) \
		$(gnomad_pop_bucket_path) \
		--freq-threshold 0.001

.PHONY: download-gnomad-out
download-gnomad-out:
	@echo "Downloading gnomAD allele frequencies..."
	@mkdir -p data/gnomad
	gsutil -m cp -r $(gnomad_af_bucket_path) data/gnomad/
	gsutil -m cp -r $(gnomad_pop_bucket_path) data/gnomad/

.PHONY: analyze-freq-dist
analyze-freq-dist:
	@echo "Analyzing frequency distribution..."
	uv run scripts/compute_variation_ratios.py \
		--vcfs-path "./data/vcf/*.vcf.bgz" \
		--gnomad-va-file ./data/gnomad/$(gnomad_af_base) \
		--gnomad-sa-file ./data/gnomad/$(gnomad_pop_base) \
		--output-ht ./data/analysis/freq_ht_per_sample.ht

.PHONY: run-haplotype-computation
run-haplotype-computation:
	@echo "Running haplotype computation..."
	@mkdir -p data/haplotypes/
	uv run scripts/compute_haplotypes.py \
		--vcfs-path "./data/vcf/*.vcf.bgz" \
		--gnomad-va-file ./data/gnomad/$(gnomad_af_base) \
		--gnomad-sa-file ./data/gnomad/$(gnomad_pop_base) \
		--window-size 100 \
		--freq-threshold 0.005 \
		--output-base ./data/haplotypes/hgdp_gnomad_merge

.PHONY: download-reference-fasta
download-reference-fasta:
	@echo "Downloading reference FASTA..."
	@mkdir -p data/reference
	gsutil -m cp gs://hail-common/references/Homo_sapiens_assembly38.fasta"*" ./data/reference

.PHONY: generate-divref
generate-divref:
	@echo "Generating divref..."
	@rm -r ./data/divref/
	@mkdir -p data/divref
	uv run scripts/create_fasta_and_index.py \
		--haplotypes-table-path ./data/haplotypes/hgdp_gnomad_merge.ht \
		--gnomad-va-file ./data/gnomad/$(gnomad_af_base) \
		--reference-fasta ./data/reference/Homo_sapiens_assembly38.fasta.gz \
		--window-size 25 \
		--output-base ./data/divref/DivRef
	rm ./data/divref/.*.crc

.PHONY: generate-merged-divref
generate-merged-divref:
	@echo "Generating divref..."
	@rm -r ./data/divref-merged/
	@mkdir -p data/divref-merged
	# ensure 32G memory for shuffles
	uv run scripts/create_fasta_and_index.py \
		--haplotypes-table-path ./data/haplotypes/hgdp_gnomad_merge.ht \
		--gnomad-va-file ./data/gnomad/$(gnomad_af_base) \
		--reference-fasta ./data/reference/Homo_sapiens_assembly38.fasta.gz \
		--window-size 25 \
		--output-base ./data/divref-merged/DivRef \
		--merge --split-contigs
	rm ./data/divref-merged/.*.crc

.PHONY: index-fasta
index-fasta:
	$(samtools) dict \
		./data/divref/DivRef.haplotypes.fasta \
		-o ./data/divref/DivRef.haplotypes.fasta.dict \
		-u '.'
	$(samtools) faidx \
		./data/divref/DivRef.haplotypes.fasta

.PHONY: index-fastas
index-fastas:
	sh ./scripts/index_all.sh ./data/divref $(samtools)


.PHONY: index-merged-fasta
index-merged-fasta:
	$(samtools) dict \
		./data/divref/DivRef.haplotypes_gnomad_merge.fasta \
		-o ./data/divref/DivRef.haplotypes_gnomad_merge.fasta.dict \
		-u '.'
	$(samtools) faidx \
		./data/divref/DivRef.haplotypes_gnomad_merge.fasta


.PHONY: create-gnomad-sites-vcf
create-gnomad-sites-vcf:
	uv run scripts/create_gnomad_sites_vcf.py \
		./data/gnomad/$(gnomad_af_base) \
		./data/gnomad/gnomad_variants.common.vcf.bgz \
		--min-popmax 0.005

.PHONY: bundle
bundle:
	rm -rf ./dist/
	mkdir -p dist/staging
	cp README.md dist/staging/
	cp LICENSE dist/staging/
	cp Citation.cff dist/staging
	cp bundle/* dist/staging
	cp -r ./data/divref-merged dist/staging/DivRef
	cp scripts/remap_divref.py dist/staging/DivRef
	cp -r ./data/divref dist/staging/DivRef_just_haplotypes
	cp scripts/remap_divref.py dist/staging/DivRef_just_haplotypes

.PHONY: run-frequency-calc
run-frequency-calc:
	uv run scripts/compute_haplotype_statistics.py \
		--haplotypes-table-path ./data/haplotypes/hgdp_gnomad_merge.ht \
		--gnomad-va-file ./data/gnomad/$(gnomad_af_base) \
		--output-base ./data/haplotypes/stats \
		--frequency-cutoffs 0.08,0.04,0.02,0.01,0.005 \
		--window-sizes 5,10,15,20,30,100


