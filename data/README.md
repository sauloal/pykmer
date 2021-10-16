# Data

## Process

```bash
K=15
for F in data/*.fa.bgz; do
    echo $F;
    G=$F.$K.kin;
    Z=$G.bgz;
    
    if [[ -f "$Z" ]]; then
        echo "$F exists";
    else
        echo "$F processing  - $K";
        pypy ./indexer.py $F $K;
        if [[ ! -f "$G" ]]; then
            echo "$F compressing - $K - $G";
            bgzip -i -I $G.bgz.gzi -l 9 -c $G > $G.bgz.tmp && mv $G.bgz.tmp $Z;
        fi
    fi
done
```

## Download

```bash
function download_fasta {
    OUT=$1
    URL=$2
    if [[ -f "${OUT}" ]]; then
        echo "${OUT} already exists"
    else
        echo "${OUT} downloading"
        wget "$URL" -O - | bgzip -c > "${OUT}.tmp" && mv "${OUT}.tmp" "${OUT}"
    fi
}

function download_fasta_gz {
    OUT=$1
    URL=$2
    if [[ -f "${OUT}" ]]; then
        echo "${OUT} already exists"
    else
        echo "${OUT} downloading"
        wget "$URL" -O - | gunzip -c | bgzip -c > "${OUT}.tmp" && mv "${OUT}.tmp" "${OUT}"
    fi
}

function download_fasta_zip {
    OUT=$1
    URL=$2
    BN=`basename ${URL}`
    FA=${BN//.zip/}
    if [[ -f "${OUT}" ]]; then
        echo "${OUT} already exists"
    else
        echo "${OUT} downloading BN ${BN} FA ${FA}"
        if [[ -f "${BN}" ]]; then rm -v ${BN}; fi
        if [[ -f "${FA}" ]]; then rm -v ${FA}; fi
        wget "$URL" -O "${BN}"
        unzip "${BN}"
        bgzip -c "${FA}" > "${OUT}.tmp" && mv "${OUT}.tmp" "${OUT}"
        rename 's/.fas.bgz/.fa.bgz/'   *.fas.bgz
        rename 's/.fasta.bgz/.fa.bgz/' *.fasta.bgz
        if [[ -f "${BN}" ]]; then rm -v ${BN}; fi
        if [[ -f "${FA}" ]]; then rm -v ${FA}; fi
    fi
}



CORNELL=ftp://ftp.sgn.cornell.edu/genomes

download_fasta     Arabidopsis_thaliana_10.fa.bgz                                          ${CORNELL}/Arabidopsis_thaliana/assembly/TAIR10/TAIR10_genome.fas
download_fasta_gz  Capsicum_annuum_UCD10X_v1.0.fa.bgz                                      ${CORNELL}/Capsicum_annuum/C.annuum_UCD10X/Capsicum_annuum_UCD10X_v1.0.fasta.gz
download_fasta_gz  Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0.fa.bgz          ${CORNELL}/Capsicum_annuum/C.annuum_glabriusculum/Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0.fasta.gz
download_fasta_gz  Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0-genes.fa.bgz    ${CORNELL}/Capsicum_annuum/C.annuum_glabriusculum/Chiltepin_genes_sequence_v2.0.fa.gz
download_fasta_gz  Capsicum.annuum.L_Zunla-1_Release_2.0.fa.bgz                            ${CORNELL}/Capsicum_annuum/C.annuum_zunla/assemblies/Capsicum.annuum.L_Zunla-1_Release_2.0.fasta.gz
download_fasta_gz  Capsicum_chinense.v0.5.fa.bgz                                           ${CORNELL}/Capsicum_chinense/assemblies/chinense.v0.5.fasta.gz
download_fasta     Coffea_humblotiana.fa.bgz                                               ${CORNELL}/Coffea_humblotiana/cohum_chrom.fasta
download_fasta_gz  Coffea_canephora.fa.bgz                                                 ${CORNELL}/Coffea_canephora/assembly/pseudomolecules.fa.gz
download_fasta     Nicotiana_attenuata_v2.fa.bgz                                           ${CORNELL}/Nicotiana_attenuata/NIATTr2.chr.fa
download_fasta_gz  Nicotiana_benthamiana_v1.0.1.scaffolds.nrcontigs.fa.bgz                 ${CORNELL}/Nicotiana_benthamiana/assemblies/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta.gz
download_fasta_gz  Nicotiana_otophora_AWOL01S.fa.bgz                                       ${CORNELL}/Nicotiana_otophora/assembly/Noto_AWOL01S.fa.gz
download_fasta_gz  Nicotiana_sylvestris_ASAF01.fa.bgz                                      ${CORNELL}/Nicotiana_sylvestris/assembly/Nsyl_ASAF01.fa.gz
download_fasta_gz  Nicotiana_tabacum_v4.5_genome_Scf_Edwards2017.fa.bgz                    ${CORNELL}/Nicotiana_tabacum/edwards_et_al_2017/assembly/Nitab-v4.5_genome_Scf_Edwards2017.fasta.gz
download_fasta_gz  Nicotiana_tomentosiformis.fa.bgz                                        ${CORNELL}/Nicotiana_tomentosiformis/assembly/Ntom_ASAG01.fa.gz
download_fasta     Oryza_sativa_v5.00_masked.fa.bgz                                        ${CORNELL}/Oryza_sativa/assembly/build_5.00/IRGSPb5.fa.masked
download_fasta     Petunia_axillaris_v1.6.2_genome.fa.bgz                                  ${CORNELL}/Petunia_axillaris/assembly/Petunia_axillaris_v1.6.2_genome.fasta
download_fasta     Petunia_inflata_v1.0.1_genome.fa.bgz                                    ${CORNELL}/Petunia_inflata/assembly/Petunia_inflata_v1.0.1_genome.fasta
download_fasta     Physalis_pruinosa-genome_v1.fa.bgz                                      ${CORNELL}/Physalis_pruinosa/assembly/build_1.00/Physalis-pruinosa-genome_v1.fa
download_fasta     Solanum_chilense.scaffolds.fa.bgz                                       ${CORNELL}/Solanum_chilense/Stam_et_al_2019/Solanum_chilense.scaffolds.fa
download_fasta     Solanum_galapagense_LA0436.fa.bgz                                       ${CORNELL}/Solanum_galapagense/LA0436/assemblies/denovo/Sg_LA0436_denovo.fa
download_fasta     Solanum_pennellii.fa.bgz                                                ${CORNELL}/Solanum_pennellii/Spenn.fasta
download_fasta     Solanum_pennellii_LYC1722_canu.fa.bgz                                   ${CORNELL}/Solanum_pennellii_LYC1722/LYC1722/canu-smartdenovo_pass5.fasta
download_fasta     Solanum_lycopersicoides_LA2951_v1.0_chromosomes.fa.bgz                  ${CORNELL}/Solanum_lycopersicoides/SlydLA2951_v1.0_chromosomes_contigs.fasta
download_fasta     Solanum_lycopersicum_cerasiforme_LA1673_r1.1.fa.bgz                     ${CORNELL}/Solanum_lycopersicum_cerasiforme/LA1673/SLYcer_r1.1.fasta
download_fasta_gz  Solanum_lycopersicum_FLA1.3.fa.bgz                                      ${CORNELL}/Solanum_lycopersicum/Fla.8924/FLA1.3_genome.fasta.gz
download_fasta_gz  Solanum_lycopersicum_M82_1.3_genome.fa.bgz                              ${CORNELL}/Solanum_lycopersicum/M82/M821.3_genome.fasta.gz
download_fasta_gz  Solanum_lycopersicum_scaffolds.1.00.fa.bgz                              ${CORNELL}/Solanum_lycopersicum/Heinz1706/wgs/assembly/build_1.00/S_lycopersicum_scaffolds.1.00.fa.gz
download_fasta_gz  Solanum_lycopersicum_chromosomes.4.00.fa.bgz                            ${CORNELL}/Solanum_lycopersicum/Heinz1706/wgs/assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa.gz
download_fasta     Solanum_melongena_HQ-1315_HQ-reference.fa.bgz                           ${CORNELL}/Solanum_melongena_HQ-1315/01.SME-HQ-reference.fasta
download_fasta_gz  Solanum_melongena_r2.5.1.fa.bgz                                         ${CORNELL}/Solanum_melongena/SME_r2.5.1.fa.gz
download_fasta     Solanum_pimpinellifolium_LA1589_2020_v0.1_chromosomes.fa.bgz            ${CORNELL}/Solanum_pimpinellifolium/LA1589/2020/SpimLA1589_v0.1_chromosomes.fa
download_fasta     Solanum_pimpinellifolium_LA1670_SPI_r1.1.pmol.fa.bgz                    ${CORNELL}/Solanum_pimpinellifolium/LA1670/SPI_r1.1.pmol.fasta
download_fasta     Solanum_pimpinellifolium_LA2093_v1.5.fa.bgz                             ${CORNELL}/Solanum_pimpinellifolium/LA2093/Spimp_LA2093_genome_v1.5/LA2093_genome_v1.5.fa 
download_fasta_gz  Solanum_pimpinellifolium_LA0480_scaffolds.abyss.77.fa.bgz               ${CORNELL}/Solanum_pimpinellifolium/LA0480/assembly/tomato-scaffolds.abyss.77.fasta.gz
download_fasta_gz  Solanum_pimpinellifolium_LA1589_1.0_contigs.fa.bgz                      ${CORNELL}/Solanum_pimpinellifolium/LA1589/assemblies/A-1.0/Spimpinellifolium_genome.contigs.fasta.gz
download_fasta_gz  Solanum_pimpinellifolium_BGV006775_BGV1.3_genome.fa.bgz                 ${CORNELL}/Solanum_pimpinellifolium/BGV006775/BGV1.3_genome.fasta.gz
download_fasta     Solanum_tuberosum_PGSC_DM_v3_superscaffolds.fa.bgz                      ${CORNELL}/Solanum_tuberosum/assembly/PGSC_DM_v3/PGSC_DM_v3_superscaffolds.fasta
download_fasta_zip Solanum_tuberosum_PGSC_DM_v4.03_pseudomolecules.fa.bgz                  ${CORNELL}/Solanum_tuberosum/assembly/PGSC_DM_v4.03/PGSC_DM_v4.03_pseudomolecules.fasta.zip
download_fasta     Vitis_vinifera_Genoscope_12X_2010_02_12_chr.fa.bgz                      ${CORNELL}/Vitis_vinifera/assembly/Genoscope_12X_2010_02_12/chr.fa


sha256sum *.fa.bgz | tee checksum.sha256sum
```
