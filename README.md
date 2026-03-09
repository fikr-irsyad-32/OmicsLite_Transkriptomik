### Capstone Project

**Judul:** Gene Expression Data from Human Peripheral Blood

**Contributor:** Chao S, Pilcz T, Stamatiou D, Ying J, Dempsey AA

**GEO accession:** GSE203024

**Pendahuluan:**

Diagnosis dini pada kanker telah terbukti meningkatkan angka harapan hidup untuk banyak tipe kanker. Dengan metodologi saat ini, diagnosis dini sulit dilakukan untuk kanker jaringan dalam. Namun, whole peripheral blood telah terbukti menjadi jaringan pengganti non-invasif yang menjanjikan untuk mendeteksi banyak jenis kanker.

**Metode:**

DNA microarray platform:

Affymetrix Human Genome U133 Plus 2.0 GeneChips

Informasi sampel:

2.485 whole peripheral blood dari relawan dengan berbagai latar belakang usia, jenis kelamin, dan status kanker

Pembagian kelompok: 

1.045 relawan dengan 12 tipe kanker yang berbeda (kanker bladder, breast, colorectal, cervical, endometrial, liver, nasopharyngeal, ovarian, prostate, stomach, kidney, dan testicular) serta satu kondisi pra-kanker (colon polyps) vs 1.800 relawan non-kanker sebagai kontrol, termasuk penderita autoimun dan kardiovaskular.

Kriteria signifikan:

Adj.P.Val (Adjusted p-value) ≤ 0,01

Tools: 

https://www.ncbi.nlm.nih.gov/geo/geo2r

R version 4.5.2

https://biit.cs.ut.ee/gprofiler/

https://www.genome.jp/kegg/mapper 

**Hasil dan Interpretasi:**

Density plot menunjukkan bahwa kelompok Cancer memiliki pola distribusi yang mirip dan tumpang tindih dengan kelompok Control (Gambar 1). Plot UMAP juga menunjukkan bahwa secara umum tidak terdapat pemisahan yang jelas antara kedua kelompok tersebut, meskipun sebagian kecil sampel Control tampak sebagai outlier (Gambar 2). Secara keseluruhan, kedua plot tersebut menunjukkan bahwa profil ekspresi global kedua kelompok relatif sebanding dan tidak memperlihatkan pola pemisahan yang mengarah pada perbedaan teknis yang mencolok. Heatmap 50 DEG teratas (berdasarkan adj.P.Val) menampilkan pola ekspresi gen pada seluruh sampel dan juga menunjukkan bahwa pemisahan klaster antara Control dan Cancer belum tampak tegas (Gambar 3).

![Density](Week%204-%20Capstone%20Project/2.%20DEG%20R/Density.png)
Gambar 1. Density Plot distribusi nilai ekspresi (log₂) seluruh gen

![UMAP](Week%204-%20Capstone%20Project/2.%20DEG%20R/UMAP.png)
Gambar 2. Plot UMAP: Control vs Cancer

![Heatmap](Week%204-%20Capstone%20Project/2.%20DEG%20R/Heatmap.png)
Gambar 3. Heatmap Top 50 DEG (berdasarkan adj.P.Val)

Volcano plot hasil analisis DEG pada gabungan 13 tipe kanker/pra-kanker dibandingkan dengan kontrol menunjukkan adanya gen yang teridentifikasi sebagai downregulated dan upregulated berdasarkan kriteria adj.P.Val ≤ 0,01 dan log₂FC ≥ 0,585 atau ≤ -0,585. Nilai log₂FC 0,585 setara dengan perubahan ekspresi gen sebesar 1,5 kali. Sebagian besar gen menunjukkan nilai adj.P.Val yang sangat rendah, yang menandakan signifikansi statistik yang tinggi. Pada ambang log₂FC 0,585, ditemukan 2 gen downregulated, yaitu GRAPL dan GABPB1-AS1, serta 3 gen upregulated, yaitu MYL9, ITGA2B, dan UTS2, dengan UTS2 muncul pada dua probe set, yaitu 220784_s_at dan 220785_at (Gambar 4). Kelima gen tersebut dapat dipertimbangkan sebagai kandidat gen untuk analisis lanjutan dalam mengevaluasi potensi relevansinya sebagai biomarker pada whole peripheral blood non-invasif.

![Volcano](Week%204-%20Capstone%20Project/2.%20DEG%20R/Volcano.png)
Gambar 4. Volcano plot DEG gabungan 13 tipe kanker/pra-kanker vs kontrol. Garis putus-putus horizontal = adj.P.Val 0,1; garis vertikal = log₂FC -0,322 dan 0,322; titik berwarna dengan label gen = log₂FC ≤ -0,585 dan ≥ 0,585

Analisis GO (Gene Ontology) pada gen-gen upregulated (Gambar 5) dan downregulated (Gambar 6) dengan kriteria adj.P.Val ≤ 0,1 dan log₂FC ≥ 0,322 atau ≤ -0,322, dengan nilai log₂FC 0,322 setara dengan perubahan ekspresi gen sebesar 1,25 kali, tidak menunjukkan term GO yang spesifik terkait kanker. Sebaliknya, analisis KEGG (Kyoto Encyclopedia of Genes and Genomes) mengidentifikasi satu pathway yang secara spesifik dikategorikan sebagai pathway kanker, yaitu Pathways in cancer (Gambar 7), dengan keterlibatan sejumlah gen upregulated dan downregulated.

![DOWN](Week%204-%20Capstone%20Project/3.%20GO%20and%20KEGG/DOWN.png)
Gambar 5. Gene Ontology pada gen downregulated

![UP](Week%204-%20Capstone%20Project/3.%20GO%20and%20KEGG/UP.png)
Gambar 6. Gene Ontology pada gen upregulated

![Pathways](Week%204-%20Capstone%20Project/3.%20GO%20and%20KEGG/Pathways%20in%20cancer.png)
Gambar 7. KEGG ‘Pathways in cancer’ (gen yang terlibat ditandai warna biru = downregulated; merah = upregulated)

**Kesimpulan:**

Secara keseluruhan, visualisasi global melalui density plot, UMAP, dan heatmap menunjukkan bahwa profil ekspresi antara gabungan 13 tipe kanker/pra-kanker dan kelompok control relatif sebanding serta tidak memperlihatkan pemisahan klaster yang tegas. Namun, analisis DEG berhasil mengidentifikasi gen yang berbeda signifikan, yaitu 2 gen downregulated (GRAPL dan GABPB1-AS1) serta 3 gen upregulated (MYL9, ITGA2B, dan UTS2). Analisis GO pada gen-gen upregulated dan downregulated tidak menunjukkan term yang spesifik terkait kanker, sedangkan analisis KEGG mengidentifikasi satu pathway yang secara spesifik dikategorikan sebagai pathway kanker, yaitu Pathways in cancer.

