#!/bin/bash
set -e  # Herhangi bir adımda hata olursa pipeline'ı güvenli şekilde durdurur

# ==============================================================================
# KLEBSIELLA PNEUMONIAE WGS BİYOINFORMATİK ANALİZ PIPELINE'I
# Kullanım: bash run_pipeline.sh
# Gereksinim: Tüm izolat (*.fasta) ve referans genom (reference.fna) bu dizinde olmalı
# ==============================================================================

REF_GENOME="reference.fna"

# Referans genom ve izolat dosyaları kontrolü
if [ ! -f "$REF_GENOME" ]; then
    echo "HATA: $REF_GENOME bulunamadı! Lütfen referans genomu bu dizine ekleyin."
    exit 1
fi

if ! ls *.fasta 1> /dev/null 2>&1; then
    echo "HATA: Bu dizinde hiç .fasta dosyası bulunamadı!"
    exit 1
fi

# --- ADIM 1: QUAST Kalite Kontrolü ---
echo ">>> ADIM 1: QUAST analizi başlatılıyor..."

quast.py *.fasta -R "$REF_GENOME" -o .

echo ">>> ADIM 1 (QUAST) tamamlandı."
echo "------------------------------------------------------------------"

# --- ADIM 2: CheckM Genom Tamamlanma ve Kontaminasyon Analizi ---
echo ">>> ADIM 2: CheckM analizi başlatılıyor..."

checkm lineage_wf -x fasta . .
checkm qa lineage.ms . -f checkm_summary.txt --tab_table

echo ">>> ADIM 2 (CheckM) tamamlandı! Özet tablo 'checkm_summary.txt' olarak kaydedildi."
echo "------------------------------------------------------------------"

# --- ADIM 3: BWA-MEM Hizalama ve Samtools Sıralama ---
echo ">>> ADIM 3: BWA hizalama işlemleri başlatılıyor..."

if [ ! -f "${REF_GENOME}.bwt" ]; then
    echo ">>> Referans genom BWA için indeksleniyor..."
    bwa index "$REF_GENOME"
fi

for fasta_dosyasi in *.fasta; do
    [ -e "$fasta_dosyasi" ] || continue
    bakteri_id=$(basename "$fasta_dosyasi" .fasta)
    
    echo -e "\n----------------------------------------"
    echo ">>> $bakteri_id BWA ile hizalanıyor..."
    
    bwa mem "$REF_GENOME" "$fasta_dosyasi" > "${bakteri_id}.sam"
    
    echo ">>> $bakteri_id SAM'dan BAM'a çevriliyor ve sıralanıyor..."
    samtools view -bS "${bakteri_id}.sam" | samtools sort -o "${bakteri_id}_sorted.bam"
    samtools index "${bakteri_id}_sorted.bam"
    rm "${bakteri_id}.sam"
    
    echo ">>> $bakteri_id BAŞARIYLA TAMAMLANDI."
done

echo -e "\n----------------------------------------"
echo ">>> ADIM 3 (BWA) EKSİKSİZ BİTTİ!"
echo "------------------------------------------------------------------"

# --- ADIM 4: MUMmer (Nucmer, Delta-filter, Show-snps, Show-coords) ---
echo ">>> ADIM 4: MUMmer karşılaştırmalı genomik analizleri başlatılıyor..."

for fasta_dosyasi in *.fasta; do
    [ -e "$fasta_dosyasi" ] || continue
    bakteri_id=$(basename "$fasta_dosyasi" .fasta)
    
    echo ">>> $bakteri_id için MUMmer (nucmer) işlemleri başlatıldı..."
    
    nucmer -p "${bakteri_id}" "$REF_GENOME" "$fasta_dosyasi"
    delta-filter -1 "${bakteri_id}.delta" > "${bakteri_id}.filter"
    show-snps -Clr "${bakteri_id}.filter" > "${bakteri_id}.snps"
    show-coords -r "${bakteri_id}.filter" > "${bakteri_id}.coords"
    
    echo ">>> $bakteri_id bitti."
done

echo ">>> ADIM 4 (MUMmer) tamamlandı!"
echo "------------------------------------------------------------------"

# --- ADIM 5: Varyant Matrisi Oluşturma ---
echo ">>> ADIM 5: Varyant matrisi oluşturuluyor..."

python3 - << 'EOF'
import os
import pandas as pd
import glob

all_snps_files = glob.glob("*.snps")
variant_data = {}

print("Varyant dosyaları taranıyor...")

for file in all_snps_files:
    strain_name = os.path.basename(file).replace(".snps", "")
    variant_data[strain_name] = {}
    
    with open(file, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        parts = line.strip().split()
        if len(parts) > 5 and parts[0].isdigit():
            ref_pos = parts[0]
            ref_base = parts[1]
            query_base = parts[2]
            
            variant_id = f"Pos_{ref_pos}_{ref_base}>{query_base}"
            variant_data[strain_name][variant_id] = 1

print("Veriler birleştiriliyor ve matris oluşturuluyor...")
df = pd.DataFrame.from_dict(variant_data, orient='index').fillna(0).astype(int)
df = df.sort_index()

output_file = "variant_matrix.csv"
df.to_csv(output_file)

print(f"\nİşlem tamam! Matris başarıyla kaydedildi: {output_file}")
print(f"Toplam Suş Sayısı: {df.shape[0]}")
print(f"Toplam Benzersiz Varyant (Feature) Sayısı: {df.shape[1]}")
EOF

echo ">>> ADIM 5 (Varyant Matrisi) tamamlandı!"
echo "------------------------------------------------------------------"

# --- ADIM 6: Prokka Gen Anotasyonu ---
echo ">>> ADIM 6: Prokka gen anotasyonu başlatılıyor..."

for fasta_dosyasi in *.fasta; do
    [ -e "$fasta_dosyasi" ] || continue
    bakteri_id=$(basename "$fasta_dosyasi" .fasta)
    
    echo ">>> Suş ${bakteri_id} Prokka ile anote ediliyor..."
    prokka --outdir . \
           --prefix "${bakteri_id}" \
           --locustag "${bakteri_id}" \
           "$fasta_dosyasi" \
           --force
done

echo ">>> ADIM 6 (Prokka) tamamlandı!"
echo "------------------------------------------------------------------"

# --- ADIM 7: Panaroo Pan-Genome Analizi ---
echo ">>> ADIM 7: Panaroo pan-genome analizi 'strict' modda başlatılıyor..."

panaroo -i *.gff \
        -o . \
        --mode strict

[ -f "gene_presence_absence.csv" ] && mv gene_presence_absence.csv panaroo_gene_presence_absence.csv
[ -f "gene_presence_absence.Rtab" ] && mv gene_presence_absence.Rtab panaroo_gene_presence_absence.Rtab

echo ">>> ADIM 7 (Panaroo) tamamlandı!"
echo "------------------------------------------------------------------"

# --- ADIM 8: Roary Pan-Genome Analizi ---
echo ">>> ADIM 8: Roary pan-genome analizi başlatılıyor..."

roary -f . -i 95 -v *.gff

[ -f "gene_presence_absence.csv" ] && mv gene_presence_absence.csv roary_gene_presence_absence.csv
[ -f "gene_presence_absence.Rtab" ] && mv gene_presence_absence.Rtab roary_gene_presence_absence.Rtab

echo ">>> ADIM 8 (Roary) tamamlandı!"
echo "------------------------------------------------------------------"

# --- ADIM 9: NCBI AMRFinderPlus Antimikrobiyal Direnç Taraması ---
echo ">>> ADIM 9: AMRFinderPlus analizi başlatılıyor..."

for fasta_dosyasi in *.fasta; do
    [ -e "$fasta_dosyasi" ] || continue
    bakteri_id=$(basename "$fasta_dosyasi" .fasta)
    
    echo ">>> $bakteri_id AMRFinderPlus ile taranıyor..."
    amrfinder -n "$fasta_dosyasi" \
              -p "${bakteri_id}.faa" \
              -g "${bakteri_id}.gff" \
              -O Klebsiella_pneumoniae \
              --plus \
              -o "${bakteri_id}_amrfinder.tsv"
done

python3 - << 'EOF'
import os
import glob
import pandas as pd

amr_files = glob.glob("*_amrfinder.tsv")
amr_data = {}

for file in amr_files:
    strain_name = os.path.basename(file).replace("_amrfinder.tsv", "")
    amr_data[strain_name] = {}
    df_tmp = pd.read_csv(file, sep="\t")
    gene_col = "Gene symbol" if "Gene symbol" in df_tmp.columns else "Element symbol"
    if gene_col in df_tmp.columns:
        for gene in df_tmp[gene_col].dropna().unique():
            amr_data[strain_name][f"AMR_{gene}"] = 1

df_amr = pd.DataFrame.from_dict(amr_data, orient="index").fillna(0).astype(int).sort_index()
df_amr.to_csv("amrfinder_matrix.csv")
print(f"AMRFinderPlus matrisi oluşturuldu: amrfinder_matrix.csv ({df_amr.shape[0]} suş, {df_amr.shape[1]} direnç geni/mutasyonu)")
EOF

echo ">>> ADIM 9 (AMRFinderPlus) tamamlandı!"
echo "------------------------------------------------------------------"

# --- ADIM 10: Abricate (7 Veritabanı: AMR, Virülans ve Plazmid Taraması) ---
echo ">>> ADIM 10: Abricate ile 7 veritabanı taranıyor..."

for db in ncbi card resfinder megares argannot vfdb plasmidfinder; do
    echo ">>> Abricate veritabanı taranıyor: $db"
    abricate --db "$db" --minid 90 --mincov 80 *.fasta > "all_strains_${db}.tab"
    abricate --summary "all_strains_${db}.tab" > "summary_${db}.tab"
done

echo ">>> ADIM 10 (Abricate) tamamlandı!"
echo "------------------------------------------------------------------"

# --- ADIM 11: Kleborate Virülans, Kapsül (K/O Locus) ve Direnç Analizi ---
echo ">>> ADIM 11: Kleborate analizi başlatılıyor..."

kleborate -a *.fasta --all -o Kleborate_Results.tsv

echo ">>> ADIM 11 (Kleborate) tamamlandı! Sonuçlar 'Kleborate_Results.tsv' dosyasına kaydedildi."
echo "------------------------------------------------------------------"

# --- ADIM 12: MLST Tip Tayini ---
echo ">>> ADIM 12: MLST sekans tip (ST) tayini başlatılıyor..."

mlst *.fasta > mlst_results.tsv

echo ">>> ADIM 12 (MLST) tamamlandı! Sonuçlar 'mlst_results.tsv' dosyasına kaydedildi."
echo "------------------------------------------------------------------"

# --- ADIM 13: Makine Öğrenmesi İçin Master Feature Matrisi Birleştirme ---
echo ">>> ADIM 13: Feature matrisleri (AMRFinderPlus + SNP + Panaroo + MLST) birleştiriliyor..."

python3 - << 'EOF'
import os
import pandas as pd

# 1. AMRFinderPlus Direnç Matrisi
df_amr = pd.read_csv("amrfinder_matrix.csv", index_col=0)
df_amr.index = df_amr.index.astype(str)

# 2. MUMmer SNP Varyant Matrisi
df_snp = pd.read_csv("variant_matrix.csv", index_col=0)
df_snp.index = df_snp.index.astype(str)

# 3. Panaroo Pan-Genome Gen Var/Yok Matrisi (Rtab formatında satırlar gen, sütunlar suştur -> transpoze alınır)
df_pan = pd.read_csv("panaroo_gene_presence_absence.Rtab", sep="\t", index_col=0).T
df_pan.index = df_pan.index.astype(str)
df_pan.columns = [f"PAN_{col}" for col in df_pan.columns]

# 4. MLST Tip Tayini (ST sütununu One-Hot Encoding ile 0/1 özniteliklerine çevir)
df_mlst_raw = pd.read_csv("mlst_results.tsv", sep="\t", header=None)
df_mlst_raw["Strain"] = df_mlst_raw[0].apply(lambda x: os.path.basename(str(x)).replace(".fasta", ""))
df_mlst_raw["ST"] = df_mlst_raw[2].astype(str)
df_mlst = pd.get_dummies(df_mlst_raw.set_index("Strain")["ST"], prefix="MLST_ST").astype(int)

# Tüm tabloları suş ID'si üzerinden birleştir ve eksik değerleri 0 yap
master_df = df_amr.join([df_snp, df_pan, df_mlst], how="outer").fillna(0).astype(int)
master_df.index.name = "Strain_ID"
master_df = master_df.sort_index()

master_df.to_csv("master_feature_matrix.csv")

print(f"\n>>> Master Feature Matrisi oluşturuldu: master_feature_matrix.csv")
print(f"Toplam Suş Sayısı: {master_df.shape[0]}")
print(f"Toplam Birleştirilmiş Öznitelik (Feature) Sayısı: {master_df.shape[1]}")
EOF

echo ">>> TÜM BİYOINFORMATİK PIPELINE BAŞARIYLA TAMAMLANDI! <<<"
