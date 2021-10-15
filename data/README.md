# Data

```bash
wget ftp://ftp.sgn.cornell.edu/genomes/Arabidopsis_thaliana/assembly/TAIR10/TAIR10_genome.fas -O - | bgzip -c > Arabidopsis_thaliana_10.fas.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Capsicum_annuum/C.annuum_UCD10X/Capsicum_annuum_UCD10X_v1.0.fasta.gz -O - | gunzip -c | bgzip -c > Capsicum_annuum_UCD10X_v1.0.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Capsicum_annuum/C.annuum_glabriusculum/Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0.fasta.gz -O - | gunzip -c | bgzip -c > Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Capsicum_annuum/C.annuum_glabriusculum/Chiltepin_genes_sequence_v2.0.fa.gz -O - | gunzip -c | bgzip -c > Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0-genes.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Capsicum_annuum/C.annuum_zunla/assemblies/Capsicum.annuum.L_Zunla-1_Release_2.0.fasta.gz -O - | gunzip -c | bgzip -c > Capsicum.annuum.L_Zunla-1_Release_2.0.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Capsicum_chinense/assemblies/chinense.v0.5.fasta.gz -O - | gunzip -c | bgzip -c > Capsicum_chinense.v0.5.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Coffea_canephora/assembly/pseudomolecules.fa.gz -O - | gunzip -c | bgzip -c > Coffea_canephora.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Coffea_humblotiana/cohum_chrom.fasta -O - | bgzip -c > Coffea_humblotiana.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Nicotiana_attenuata/NIATTr2.chr.fa -O Nicotiana_attenuata_v2.fasta -O - | bgzip -c > Nicotiana_attenuata_v2.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Nicotiana_benthamiana/assemblies/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta.gz -O - | gunzip -c | bgzip -c > Nicotiana_benthamiana_v1.0.1.scaffolds.nrcontigs.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Nicotiana_otophora/assembly/Noto_AWOL01S.fa.gz -O - | gunzip -c | bgzip -c > Nicotiana_otophora_AWOL01S.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Nicotiana_sylvestris/assembly/Nsyl_ASAF01.fa.gz -O - | gunzip -c | bgzip -c > Nicotiana_sylvestris_ASAF01.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Nicotiana_tabacum/edwards_et_al_2017/assembly/Nitab-v4.5_genome_Scf_Edwards2017.fasta.gz -O - | gunzip -c | bgzip -c > Nicotiana_tabacum_v4.5_genome_Scf_Edwards2017.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Nicotiana_tomentosiformis/assembly/Ntom_ASAG01.fa.gz -O - | gunzip -c | bgzip -c > Nicotiana_tomentosiformis.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Oryza_sativa/assembly/build_5.00/IRGSPb5.fa.masked -O Oryza_sativa_v5.00_masked.fa -O - | bgzip -c > Oryza_sativa_v5.00_masked.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Petunia_axillaris/assembly/Petunia_axillaris_v1.6.2_genome.fasta -O - | bgzip -c > Petunia_axillaris_v1.6.2_genome.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Petunia_inflata/assembly/Petunia_inflata_v1.0.1_genome.fasta -O - | bgzip -c > Petunia_inflata_v1.0.1_genome.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Physalis_pruinosa/assembly/build_1.00/Physalis-pruinosa-genome_v1.fa -O - | bgzip -c > Physalis_pruinosa-genome_v1.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_chilense/Stam_et_al_2019/Solanum_chilense.scaffolds.fa -O - | bgzip -c > Solanum_chilense.scaffolds.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_galapagense/LA0436/assemblies/denovo/Sg_LA0436_denovo.fa -O - | bgzip -c > Solanum_galapagense_LA0436.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_lycopersicoides/SlydLA2951_v1.0_chromosomes_contigs.fasta -O - | bgzip -c > Solanum_lycopersicoides_LA2951_v1.0_chromosomes_contigs.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_lycopersicum/Fla.8924/FLA1.3_genome.fasta.gz -O - | gunzip -c | bgzip -c > Solanum_lycopersicum_FLA1.3.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_lycopersicum/M82/M821.3_genome.fasta.gz -O - | gunzip -c | bgzip -c > Solanum_lycopersicum_M82_1.3_genome.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_lycopersicum_cerasiforme/LA1673/SLYcer_r1.1.fasta -O - | bgzip -c > Solanum_lycopersicum_cerasiforme_LA1673_r1.1.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_melongena/SME_r2.5.1.fa.gz -O - | gunzip -c | bgzip -c > Solanum_melongena_r2.5.1.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_melongena_HQ-1315/01.SME-HQ-reference.fasta -O - | bgzip -c > Solanum_melongena_HQ-1315_HQ-reference.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pennellii/Spenn.fasta -O - | bgzip -c > Solanum_pennellii.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pennellii_LYC1722/LYC1722/canu-smartdenovo_pass5.fasta -O - | bgzip -c > Solanum_pennellii_LYC1722_canu-smartdenovo_pass5.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_peruvianum/Speru_denovo.fa -O - | bgzip -c > Solanum_peruvianum.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_tuberosum/assembly/PGSC_DM_v3/PGSC_DM_v3_superscaffolds.fasta -O - | bgzip -c > Solanum_tuberosum_PGSC_DM_v3_superscaffolds.fasta.bgzip

wget ftp://ftp.sgn.cornell.edu/genomes/Vitis_vinifera/assembly/Genoscope_12X_2010_02_12/chr.fa -O - | bgzip -c > Vitis_vinifera_Genoscope_12X_2010_02_12_chr.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pimpinellifolium/LA0480/assembly/tomato-scaffolds.abyss.77.fasta.gz -O - | gunzip -c | bgzip -c > Solanum_pimpinellifolium_LA0480_scaffolds.abyss.77.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pimpinellifolium/LA1589/assemblies/A-1.0/Spimpinellifolium_genome.contigs.fasta.gz -O - | gunzip -c | bgzip -c > Solanum_pimpinellifolium_LA1589_1.0_contigs.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pimpinellifolium/LA1589/2020/SpimLA1589_v0.1_chromosomes.fa -O - | bgzip -c > Solanum_pimpinellifolium_LA1589_2020_v0.1_chromosomes.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pimpinellifolium/LA1670/SPI_r1.1.pmol.fasta -O - | bgzip -c > Solanum_pimpinellifolium_LA1670_SPI_r1.1.pmol.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pimpinellifolium/LA2093/Spimp_LA2093_genome_v1.5/LA2093_genome_v1.5.fa -O - | bgzip -c > Solanum_pimpinellifolium_LA2093_v1.5.fa.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_pimpinellifolium/BGV006775/BGV1.3_genome.fasta.gz -O - | gunzip -c | bgzip -c > Solanum_pimpinellifolium_BGV006775_BGV1.3_genome.fasta.bgz

wget ftp://ftp.sgn.cornell.edu/genomes/Solanum_tuberosum/assembly/PGSC_DM_v4.03/PGSC_DM_v4.03_pseudomolecules.fasta.zip 


for GZ in *.gz; do
    BGZ=${GZ//.gz/.bgz}
    echo "$GZ $BGZ"
    gunzip -c $GZ | bgzip -c > $BGZ
    rm -v $GZ
done
```
