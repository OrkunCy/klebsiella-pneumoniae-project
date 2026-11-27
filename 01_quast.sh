#!/bin/bash

# Adım 1: Tüm Genomların QUAST Kalite Kontrol Analizi
# Kullanım: bash 01_quast.sh
# Gereksinim: Örnek izolatlar (*.fasta) ve referans genom (reference.fna) bu script ile aynı dizinde olmalıdır.

REF_GENOME="reference.fna"
OUTPUT_DIR="QUAST_TUM_GENOMLAR_RAPORU"

if [ ! -f "$REF_GENOME" ]; then
    echo "HATA: Referans genom dosyası ($REF_GENOME) bulunamadı!"
    echo "Lütfen referans genomu '$REF_GENOME' adıyla bu klasöre ekleyin."
    exit 1
fi

echo "QUAST analizi başlatılıyor..."

quast.py *.fasta \
    -R "$REF_GENOME" \
    -o "$OUTPUT_DIR"

echo "QUAST tamamlandı! Sonuçlar '$OUTPUT_DIR' klasörüne kaydedildi."
