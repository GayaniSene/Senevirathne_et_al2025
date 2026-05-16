#bedtools intersect for HAR-overlapping regions 
wc -l *bed
wc -l *narrowPeak


#activate the environment
source activate singlecell2


#first step is to get each peak file and overlap it with HARS. 
#Here, -wa retains the peaks from the first .bed file when intersected with the second file. 

bedtools intersect -a ChondroProg_peaks.narrowPeak -b HARS_hg38.bed -wa > ChondroProg_commonpeaks.bed
bedtools intersect -a EndothelialCells_peaks.narrowPeak -b HARS_hg38.bed -wa > EndothelialCells_commonpeaks.bed
bedtools intersect -a Fibro1_peaks.narrowPeak -b HARS_hg38.bed -wa > Fibro1_commonpeaks.bed
bedtools intersect -a Fibro2_peaks.narrowPeak -b HARS_hg38.bed -wa > Fibro2_commonpeaks.bed
bedtools intersect -a Fibro3_peaks.narrowPeak -b HARS_hg38.bed -wa > Fibro3_commonpeaks.bed
bedtools intersect -a ImmuneCells_peaks.narrowPeak -b HARS_hg38.bed -wa > ImmuneCells_commonpeaks.bed
bedtools intersect -a Mes1_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes1_commonpeaks.bed
bedtools intersect -a Mes2_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes2_commonpeaks.bed
bedtools intersect -a MyoProg+Pax3_peaks.narrowPeak -b HARS_hg38.bed -wa > MyoProg+Pax3_commonpeaks.bed
bedtools intersect -a MyoProg_peaks.narrowPeak -b HARS_hg38.bed -wa > MyoProg_commonpeaks.bed
bedtools intersect -a Perichondro+Osteoblasts_peaks.narrowPeak -b HARS_hg38.bed -wa > Perichondro+Osteoblasts_commonpeaks.bed
bedtools intersect -a Perimysium_peaks.narrowPeak -b HARS_hg38.bed -wa > Perimysium_commonpeaks.bed
bedtools intersect -a ProxMes_peaks.narrowPeak -b HARS_hg38.bed -wa > ProxMes_commonpeaks.bed
bedtools intersect -a RestingChondro2_peaks.narrowPeak -b HARS_hg38.bed -wa > RestingChondro2_commonpeaks.bed
bedtools intersect -a RestingChondro_peaks.narrowPeak -b HARS_hg38.bed -wa > RestingChondro_commonpeaks.bed
bedtools intersect -a Schwann_peaks.narrowPeak -b HARS_hg38.bed -wa > Schwann_commonpeaks.bed
bedtools intersect -a SmoothMuscleProg_peaks.narrowPeak -b HARS_hg38.bed -wa > SmoothMuscleProg_commonpeaks.bed
bedtools intersect -a Teno_peaks.narrowPeak -b HARS_hg38.bed -wa > Teno_peaks_commonpeaks.bed
bedtools intersect -a TransMes_peaks.narrowPeak -b HARS_hg38.bed -wa > TransMes_commonpeaks.bed



#concatenate overlapping files to look at unique peaks

cat EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_chondroprog.bed

cat ChondroProg_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_endothelial.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_fibro1.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_fibro2.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_fibro3.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed  Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_immunecells.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_mes1.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_mes2.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_myoprog.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_myoprogPax7.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_perichondroOsteo.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed  ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_perimysium.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_proxmes.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total_restingchondro.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total-restingchondro2.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total-schwann.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed Teno_peaks_commonpeaks.bed TransMes_commonpeaks.bed >> total-smooth.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed TransMes_commonpeaks.bed >> total-teno.bed

cat ChondroProg_commonpeaks.bed EndothelialCells_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Fibro3_commonpeaks.bed ImmuneCells_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoProg_commonpeaks.bed MyoProg+Pax3_commonpeaks.bed Perichondro+Osteoblasts_commonpeaks.bed Perimysium_commonpeaks.bed ProxMes_commonpeaks.bed RestingChondro_commonpeaks.bed RestingChondro2_commonpeaks.bed Schwann_commonpeaks.bed SmoothMuscleProg_commonpeaks.bed Teno_peaks_commonpeaks.bed >> total-transMes.bed




#intersect file A (X cell type overlapping with HARs) and file B (total HARs overlapping peaks-X) to retain only the peaks that's unique (for X cell) in file A. (for this use #-v command)



bedtools intersect -a ChondroProg_commonpeaks.bed -b total_chondroprog.bed -v > ChondroProg_unipeaks.bed
bedtools intersect -a EndothelialCells_commonpeaks.bed -b total_endothelial.bed -v > EndothelialCells_unipeaks.bed
bedtools intersect -a Fibro1_commonpeaks.bed -b total_fibro1.bed -v > Fibro1_unipeaks.bed
bedtools intersect -a Fibro2_commonpeaks.bed -b total_fibro2.bed  -v > Fibro2_unipeaks.bed
bedtools intersect -a Fibro3_commonpeaks.bed -b total_fibro3.bed -v > Fibro3_unipeaks.bed
bedtools intersect -a ImmuneCells_commonpeaks.bed -b total_immunecells.bed -v > ImmuneCells_unipeaks.bed
bedtools intersect -a Mes1_commonpeaks.bed -b total_mes1.bed -v > Mes1_unipeaks.bed
bedtools intersect -a Mes2_commonpeaks.bed -b total_mes2.bed -v > Mes2_unipeaks.bed
bedtools intersect -a MyoProg+Pax3_commonpeaks.bed -b total_myoprogPax7.bed -v > MyoProg+Pax3_unipeaks.bed
bedtools intersect -a MyoProg_commonpeaks.bed -b total_myoprog.bed -v > MyoProg_unipeaks.bed
bedtools intersect -a Perichondro+Osteoblasts_commonpeaks.bed -b total_perichondroOsteo.bed -v > Perichondro+Osteoblasts_unipeaks.bed
bedtools intersect -a Perimysium_commonpeaks.bed -b total_perimysium.bed -v > Perimysium_unipeaks.bed
bedtools intersect -a ProxMes_commonpeaks.bed -b total_proxmes.bed -v > ProxMes_unipeaks.bed
bedtools intersect -a RestingChondro2_commonpeaks.bed -b total-restingchondro2.bed -v > RestingChondro2_unipeaks.bed
bedtools intersect -a RestingChondro_commonpeaks.bed -b total_restingchondro.bed -v > RestingChondro_unipeaks.bed
bedtools intersect -a Schwann_commonpeaks.bed -b total-schwann.bed -v > Schwann_unipeaks.bed
bedtools intersect -a SmoothMuscleProg_commonpeaks.bed -b total-smooth.bed -v > SmoothMuscleProg_unipeaks.bed
bedtools intersect -a Teno_peaks_commonpeaks.bed -b total-teno.bed -v > Teno_peaks_unipeaks.bed
bedtools intersect -a TransMes_commonpeaks.bed -b total-transMes.bed -v > TransMes_unipeaks.bed



-------------------- E57 --------------


wc -l *bed


bedtools intersect -a ArterialEndo_peaks.narrowPeak -b HARS_hg38.bed -wa > ArterialEndo_commonpeaks.bed
bedtools intersect -a Erythro_peaks.narrowPeak -b HARS_hg38.bed -wa > Erythro_peaks_commonpeaks.bed
bedtools intersect -a Fibro1_peaks.narrowPeak -b HARS_hg38.bed -wa > Fibro1_commonpeaks.bed
bedtools intersect -a Fibro2_peaks.narrowPeak -b HARS_hg38.bed -wa > Fibro2_commonpeaks.bed
bedtools intersect -a Macrophages_peaks.narrowPeak -b HARS_hg38.bed -wa > Macrophages_commonpeaks.bed
bedtools intersect -a Mes1_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes1_commonpeaks.bed
bedtools intersect -a Mes2_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes2_commonpeaks.bed
bedtools intersect -a MyoC_peaks.narrowPeak -b HARS_hg38.bed -wa > MyoC_commonpeaks.bed
bedtools intersect -a perichondrium_Osteoblast_peaks.narrowPeak -b HARS_hg38.bed -wa > perichondrium_Osteoblast_commonpeaks.bed
bedtools intersect -a Perichondrium_peaks.narrowPeak -b HARS_hg38.bed -wa > Perichondrium_commonpeaks.bed
bedtools intersect -a Perimysium_peaks.narrowPeak -b HARS_hg38.bed -wa > Perimysium_commonpeaks.bed
bedtools intersect -a Myoprog+Pax7_peaks.narrowPeak -b HARS_hg38.bed -wa > Myoprog+Pax7_peaks_commonpeaks.bed
bedtools intersect -a PrimaryErythro_peaks.narrowPeak -b HARS_hg38.bed -wa > PrimaryErythro_commonpeaks.bed
bedtools intersect -a Proliferating2_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating2_commonpeaks.bed
bedtools intersect -a Proliferating3_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating3_commonpeaks.bed
bedtools intersect -a Proliferating_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating_commonpeaks.bed
bedtools intersect -a RestingChondro_peaks.narrowPeak -b HARS_hg38.bed -wa > RestingChondro_commonpeaks.bed
bedtools intersect -a Schwann_peaks.narrowPeak -b HARS_hg38.bed -wa > Schwann_commonpeaks.bed
bedtools intersect -a TransMes_peaks.narrowPeak -b HARS_hg38.bed -wa > TransMes_commonpeaks.bed
bedtools intersect -a VenousEndo_peaks.narrowPeak -b HARS_hg38.bed -wa > VenousEndo_commonpeaks.bed





cat RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed PrimaryErythro_commonpeaks.bed Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed >> total-venousendo.bed 

cat Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed PrimaryErythro_commonpeaks.bed Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed >> total-restingchondro.bed 

cat TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed PrimaryErythro_commonpeaks.bed Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed >> total-schwann.bed 

cat Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed PrimaryErythro_commonpeaks.bed Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed >> total-transMes.bed 

cat Proliferating3_commonpeaks.bed PrimaryErythro_commonpeaks.bed Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed >> total-proliferating.bed 

cat PrimaryErythro_commonpeaks.bed Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed >> total-proliferating3.bed 

cat Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed >> total-primaryerythro.bed 

cat Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed >> total-proliferating2.bed 

cat Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed >> total-MyoPax7.bed 

cat Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed >> total-Perimysium.bed 

cat MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed >> total-Mes2.bed 

cat Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed >> total-MyoC.bed 

cat perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed >> total-perichondro.bed 

cat Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed >> total-perichondroOsteo.bed 

cat Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed >> total-Macrophages.bed 

cat Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed >> total-Mes1.bed 

cat Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed >> total-Fibro1.bed 

cat Erythro_peaks_commonpeaks.bed ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed >> total-Fibro2.bed 

cat ArterialEndo_commonpeaks.bed VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed >> total-Erythro.bed 

cat VenousEndo_commonpeaks.bed RestingChondro_commonpeaks.bed Schwann_commonpeaks.bed TransMes_commonpeaks.bed Proliferating_commonpeaks.bed Proliferating3_commonpeaks.bed  PrimaryErythro_commonpeaks.bed  Proliferating2_commonpeaks.bed Myoprog+Pax7_peaks_commonpeaks.bed Perimysium_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed Perichondrium_commonpeaks.bed perichondrium_Osteoblast_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Fibro1_commonpeaks.bed Fibro2_commonpeaks.bed Erythro_peaks_commonpeaks.bed >> total-ArterialEndo.bed 


bedtools intersect -a ArterialEndo_commonpeaks.bed -b total-ArterialEndo.bed -v > ArterialEndo_unipeaks.bed
bedtools intersect -a Erythro_peaks_commonpeaks.bed -b total-Erythro.bed  -v > Erythro_peaks_unipeaks.bed
bedtools intersect -a Fibro1_commonpeaks.bed -b total-Fibro1.bed -v > Fibro1_unipeaks.bed
bedtools intersect -a Macrophages_commonpeaks.bed -b total-Macrophages.bed  -v > Macrophages_unipeaks.bed
bedtools intersect -a Mes1_commonpeaks.bed -b total-Mes1.bed  -v > Mes1_unipeaks.bed
bedtools intersect -a Mes2_commonpeaks.bed -b total-Mes2.bed  -v > Mes2_unipeaks.bed
bedtools intersect -a MyoC_commonpeaks.bed -b total-MyoC.bed  -v > MyoC_unipeaks.bed
bedtools intersect -a perichondrium_Osteoblast_commonpeaks.bed -b total-perichondroOsteo.bed  -v > perichondrium_Osteoblast_unipeaks.bed
bedtools intersect -a Perichondrium_commonpeaks.bed -b total-perichondro.bed  -v > Perichondrium_unipeaks.bed
bedtools intersect -a Perimysium_commonpeaks.bed -b total-Perimysium.bed  -v > Perimysium_unipeaks.bed
bedtools intersect -a Myoprog+Pax7_peaks_commonpeaks.bed -b total-MyoPax7.bed  -v > Myoprog+Pax7_peaks_unipeaks.bed
bedtools intersect -a PrimaryErythro_commonpeaks.bed -b total-primaryerythro.bed  -v > PrimaryErythro_unipeaks.bed
bedtools intersect -a Proliferating2_commonpeaks.bed -b total-proliferating2.bed  -v > Proliferating2_unipeaks.bed
bedtools intersect -a Proliferating3_commonpeaks.bed -b total-proliferating3.bed  -v > Proliferating3_unipeaks.bed
bedtools intersect -a Proliferating_commonpeaks.bed -b total-proliferating.bed  -v > Proliferating_unipeaks.bed
bedtools intersect -a RestingChondro_commonpeaks.bed -b total-restingchondro.bed  -v > RestingChondro_unipeaks.bed
bedtools intersect -a Schwann_commonpeaks.bed -b total-schwann.bed  -v > Schwann_unipeaks.bed
bedtools intersect -a TransMes_commonpeaks.bed -b total-transMes.bed  -v > TransMes_unipeaks.bed
bedtools intersect -a VenousEndo_commonpeaks.bed -b total-venousendo.bed  -v > VenousEndo_unipeaks.bed





----------------- E67 -------------------- 





bedtools intersect -a ArterialEndo_peaks.narrowPeak -b HARS_hg38.bed -wa > ArterialEndo_commonpeaks.bed
bedtools intersect -a Fibro1_peaks.narrowPeak -b HARS_hg38.bed -wa > Fibro1_commonpeaks.bed
bedtools intersect -a Macrophages_peaks.narrowPeak -b HARS_hg38.bed -wa > Macrophages_commonpeaks.bed
bedtools intersect -a Mes1_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes1_commonpeaks.bed
bedtools intersect -a Mes2_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes2_commonpeaks.bed
bedtools intersect -a Myocytes_peaks.narrowPeak -b HARS_hg38.bed -wa > MyoC_commonpeaks.bed
bedtools intersect -a MyoProg+7_peaks.narrowPeak -b HARS_hg38.bed -wa > MyoProg+7_peaks_commonpeaks.bed
bedtools intersect -a Perichondrium_peaks.narrowPeak -b HARS_hg38.bed -wa > Perichondrium_commonpeaks.bed
bedtools intersect -a Perimysium_peaks.narrowPeak -b HARS_hg38.bed -wa > Perimysium_commonpeaks.bed
bedtools intersect -a PrimaryErythro_peaks.narrowPeak -b HARS_hg38.bed -wa > PrimaryErythro_commonpeaks.bed
bedtools intersect -a Proliferating+Runx2_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating+Runx2_commonpeaks.bed
bedtools intersect -a Proliferating1_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating_commonpeaks.bed
bedtools intersect -a Proliferating2_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating2_commonpeaks.bed
bedtools intersect -a Proliferating3_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating3_commonpeaks.bed
bedtools intersect -a Proliferating4_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating4_commonpeaks.bed
bedtools intersect -a RestingChondro_peaks.narrowPeak -b HARS_hg38.bed -wa > RestingChondro_commonpeaks.bed
bedtools intersect -a TransMes_peaks.narrowPeak -b HARS_hg38.bed -wa > TransMes_commonpeaks.bed


wc -l *bed
cat Fibro1_commonpeaks.bed Macrophages_commonpeaks.bed Mes1_commonpeaks.bed Mes2_commonpeaks.bed MyoC_commonpeaks.bed MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_Arterial.bed

cat ArterialEndo_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_fibro1.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_macrophages.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_mes1.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_mes2.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_myoC.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_myoProg.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed  >> total_perichondrium.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_perimysium.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_primErythro.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_proliferating.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_proliferating+RUNX2.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_proliferating2.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_proliferating3.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed >> total_proliferating4.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	TransMes_commonpeaks.bed >> total_RestChondro.bed

cat ArterialEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoC_commonpeaks.bed	MyoProg+7_peaks_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed	Proliferating_commonpeaks.bed	Proliferating+Runx2_commonpeaks.bed	Proliferating2_commonpeaks.bed	Proliferating3_commonpeaks.bed	Proliferating4_commonpeaks.bed	RestingChondro_commonpeaks.bed >> total_TransMes.bed




wc -l *bed

bedtools intersect -a ArterialEndo_commonpeaks.bed -b total_Arterial.bed -v > Arterial_Endo_uni.bed	
bedtools intersect -a Fibro1_commonpeaks.bed  -b total_fibro1.bed  -v > fibro1_uni.bed	
bedtools intersect -a  Macrophages_commonpeaks.bed	 -b	 total_macrophages.bed	 -v >	 macrophages_uni.bed	
bedtools intersect -a Mes1_commonpeaks.bed	  -b total_mes1.bed   -v  >  mes1_uni.bed	
bedtools intersect -a  Mes2_commonpeaks.bed	 -b  total_mes2.bed	 -v >	mes2_uni.bed	
bedtools intersect -a  MyoC_commonpeaks.bed	 -b	 total_myoC.bed	 -v >	myoc_uni.bed	
bedtools intersect -a  MyoProg+7_peaks_commonpeaks.bed	 -b	 total_myoProg.bed	 -v >	myoprog_uni.bed	
bedtools intersect	 -a  Perichondrium_commonpeaks.bed	 -b	 total_perichondrium.bed	 -v >	perichondro_uni.bed	
bedtools intersect	 -a  Perimysium_commonpeaks.bed	 -b	 total_perimysium.bed	 -v >	perimysium_uni.bed	
bedtools intersect	 -a  PrimaryErythro_commonpeaks.bed	 -b	 total_primErythro.bed	 -v >	primErythro_uni.bed	
bedtools intersect	 -a  Proliferating2_commonpeaks.bed	 -b	 total_proliferating2.bed	 -v >	proliferating2_uni.bed	
bedtools intersect	 -a  Proliferating3_commonpeaks.bed	 -b	 total_proliferating3.bed	 -v >	proliferating3_uni.bed	
bedtools intersect	 -a  Proliferating4_commonpeaks.bed	 -b	 total_proliferating4.bed	 -v >	proliferating4_uni.bed	
bedtools intersect	 -a  Proliferating_commonpeaks.bed	 -b	 total_proliferating.bed	 -v >	proliferating_uni.bed	
bedtools intersect	 -a  Proliferating+Runx2_commonpeaks.bed -b	 total_proliferating+RUNX2.bed	 -v >	proliferating+RUNX2_uni.bed	
bedtools intersect	 -a  RestingChondro_commonpeaks.bed	 -b	 total_RestChondro.bed	 -v >	RestingChondro_uni.bed	
bedtools intersect	 -a  TransMes_commonpeaks.bed	 -b	 total_TransMes.bed	 -v >	TransMes_uni.bed	







---------------------------------- E72 ------------



E72_IL unique peaks_commonpeaks


#bedtools intersect for HAR-overlapping regions
wc -l *bed
wc -l *narrowPeak




bedtools intersect -a Erythrocytes_peaks.narrowPeak -b HARS_hg38.bed -wa > Erythrocytes_commonpeaks.bed
bedtools intersect -a Fibro1_peaks.narrowPeak -b HARS_hg38.bed -wa > Fibro1_commonpeaks.bed
bedtools intersect -a Macrophages_peaks.narrowPeak -b HARS_hg38.bed -wa > Macrophages_commonpeaks.bed
bedtools intersect -a Mes1_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes1_commonpeaks.bed
bedtools intersect -a Mes2_peaks.narrowPeak -b HARS_hg38.bed -wa > Mes2_commonpeaks.bed
bedtools intersect -a Myocytes_peaks.narrowPeak -b HARS_hg38.bed -wa > Myocytes_commonpeaks.bed
bedtools intersect -a MyoProg_peaks.narrowPeak -b HARS_hg38.bed -wa > MyoProg_commonpeaks.bed
bedtools intersect -a Myoprog+Pax7_peaks.narrowPeak -b HARS_hg38.bed -wa > Myoprog+Pax7_commonpeaks.bed
bedtools intersect -a Neurons_peaks.narrowPeak -b HARS_hg38.bed -wa > Neurons_commonpeaks.bed
bedtools intersect -a Perichondrium_peaks.narrowPeak -b HARS_hg38.bed -wa > Perichondrium_commonpeaks.bed
bedtools intersect -a Perimysium_peaks.narrowPeak -b HARS_hg38.bed -wa > Perimysium_commonpeaks.bed
bedtools intersect -a Perimysium2_peaks.narrowPeak -b HARS_hg38.bed -wa > Perimysium2_commonpeaks.bed
bedtools intersect -a PrimaryErythro_peaks.narrowPeak -b HARS_hg38.bed -wa > PrimaryErythro_commonpeaks.bed
bedtools intersect -a Proliferating1_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating1_commonpeaks.bed
bedtools intersect -a Proliferating2_peaks.narrowPeak -b HARS_hg38.bed -wa > Proliferating2_commonpeaks.bed
bedtools intersect -a RestingChondro_peaks.narrowPeak -b HARS_hg38.bed -wa > RestingChondro_commonpeaks.bed
bedtools intersect -a Schwann_peaks.narrowPeak -b HARS_hg38.bed -wa > Schwann_commonpeaks.bed
bedtools intersect -a TransMes_peaks.narrowPeak -b HARS_hg38.bed -wa > TransMes_commonpeaks.bed
bedtools intersect -a VenousEndo_peaks.narrowPeak -b HARS_hg38.bed -wa > VenousEndo_commonpeaks.bed




cat Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Proliferating1.bed	

cat Proliferating1_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Proliferating2.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-RestingChondro.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Schwann.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-TransMes.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-VenousEndo.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Erythrocytes.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Fibro1.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Macrophages.bed

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Mes1.bed

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Mes2.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Myocytes.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-MyoProg.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Myoprog+Pax7.bed

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Neurons.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Perichondrium.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium2_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Perimysium.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	PrimaryErythro_commonpeaks.bed >> total-Perimysium2.bed	

cat Proliferating1_commonpeaks.bed	Proliferating2_commonpeaks.bed	RestingChondro_commonpeaks.bed	Schwann_commonpeaks.bed	TransMes_commonpeaks.bed	VenousEndo_commonpeaks.bed	Erythrocytes_commonpeaks.bed	Fibro1_commonpeaks.bed	Macrophages_commonpeaks.bed	Mes1_commonpeaks.bed	Mes2_commonpeaks.bed	Myocytes_commonpeaks.bed	MyoProg_commonpeaks.bed	Myoprog+Pax7_commonpeaks.bed	Neurons_commonpeaks.bed	Perichondrium_commonpeaks.bed	Perimysium_commonpeaks.bed	Perimysium2_commonpeaks.bed	 >> total-PrimaryErythro.bed





bedtools intersect -a Erythrocytes_commonpeaks.bed -b total-Erythrocytes.bed -v > Erythrocytes_uni.bed
bedtools intersect -a Fibro1_commonpeaks.bed -b total-Fibro1.bed -v > Fibro1_uni.bed
bedtools intersect -a Macrophages_commonpeaks.bed -b total-Macrophages.bed -v > Macrophages_uni.bed
bedtools intersect -a Mes1_commonpeaks.bed -b total-Mes1.bed -v > uni-Mes1.bed
bedtools intersect -a Mes2_commonpeaks.bed -b total-Mes1.bed -v > uni-Mes1.bed	
bedtools intersect -a Myocytes_commonpeaks.bed -b total-Myocytes.bed -v > uni-Myocytes.bed	
bedtools intersect -a MyoProg_commonpeaks.bed -b total-MyoProg.bed -v > uni-MyoProg.bed	
bedtools intersect -a Myoprog+Pax7_commonpeaks.bed -b total-Myoprog+Pax7.bed -v > uni-Myoprog+Pax7.bed
bedtools intersect -a Neurons_commonpeaks.bed -b total-Neurons.bed -v > uni-Neurons.bed	
bedtools intersect -a Perichondrium_commonpeaks.bed -b total-Perichondrium.bed -v > uni-Perichondrium.bed	
bedtools intersect -a Perimysium2_commonpeaks.bed -b total-Perimysium2.bed -v > uni-Perimysium2.bed 	
bedtools intersect -a Perimysium_commonpeaks.bed -b total-Perimysium.bed -v > uni-Perimysium.bed	
bedtools intersect -a PrimaryErythro_commonpeaks.bed -b total-PrimaryErythro.bed -v > uni-PrimaryErythro.bed 	
bedtools intersect -a Proliferating1_commonpeaks.bed -b total-Proliferating1.bed -v > uni-Proliferating1.bed	
bedtools intersect -a Proliferating2_commonpeaks.bed -b total-Proliferating2.bed -v > uni-Proliferating2.bed
bedtools intersect -a RestingChondro_commonpeaks.bed -b total-RestingChondro.bed -v > uni-RestingChondro.bed	
bedtools intersect -a Schwann_commonpeaks.bed -b total-Schwann.bed -v > uni-Schwann.bed	
bedtools intersect -a TransMes_commonpeaks.bed -b total-TransMes.bed -v > uni-TransMes.bed	
bedtools intersect -a VenousEndo_commonpeaks.bed -b total-VenousEndo.bed -v > uni-VenousEndo.bed	





------------------------
#eRegulon motif analysis


bedtools intersect -a RestingChondro2_motifs_E53.bed -b HARS_hg38.bed -wa > RestingChondro2_motifs_HARS.bed
bedtools intersect -a Perichondro_Motif_E53.bed -b HARS_hg38.bed -wa > Perichondro_motifs_HARS.bed
bedtools intersect -a Mes2_motifs_E53.bed -b HARS_hg38.bed -wa > Mes2_motifs_HARS.bed
bedtools intersect -a Mes1_motifs_E53.bed -b HARS_hg38.bed -wa > Mes1_motifs_HARS.bed
bedtools intersect -a fibro2_motifs_E53.bed -b HARS_hg38.bed -wa > fibro2_motifs_HARS.bed
bedtools intersect -a ChondroProg_motifs_E53.bed -b HARS_hg38.bed -wa > ChondroProg_motifs_HARS.bed
bedtools intersect -a RestingChondro_motifs_E53.bed -b HARS_hg38.bed -wa > RestingChondro_motifs_HARS.bed
bedtools intersect -a Fibro1_motifs_E57.bed -b HARS_hg38.bed -wa > Fibro1_motifs_E57_HARS.bed
bedtools intersect -a RestingChondro_motifs_E57.bed -b HARS_hg38.bed -wa > RestingChondro_motifs_E57_HARS.bed
bedtools intersect -a Proliferating3_motifs_E57.bed -b HARS_hg38.bed -wa > Proliferating3_motifs_E57_HARS.bed
bedtools intersect -a Proliferating2_motifs_E57.bed -b HARS_hg38.bed -wa > Proliferating2_motifs_E57_HARS.bed
bedtools intersect -a Proliferating1_motifs_E57.bed -b HARS_hg38.bed -wa > Proliferating1_motifs_E57_HARS.bed
bedtools intersect -a Perichondrium_motifs_E57.bed -b HARS_hg38.bed -wa > Perichondrium_motifs_E57_HARS.bed
bedtools intersect -a Osteoblast_motifs_E57.bed -b HARS_hg38.bed -wa > Osteoblast_motifs_E57_HARS.bed
bedtools intersect -a Mes2_motifs_E57.bed -b HARS_hg38.bed -wa > Mes2_motifs_E57_HARS.bed






