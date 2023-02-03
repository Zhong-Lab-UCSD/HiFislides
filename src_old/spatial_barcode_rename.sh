OLD_DIR=/mnt/extraids/SDSC_NFS/linpei/hifi/recycled_flowcell/AAANK57HV
FILENAME=L1R1Dedup_AAANK57HV_21

NEW_DIR=/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANK57HV_2

mkdir -p $NEW_DIR
mkdir -p $NEW_DIR/bwa_index_L1R1

cp $OLD_DIR/$FILENAME.fasta $NEW_DIR/L1R1_dedup.fasta
cp $OLD_DIR/$FILENAME.read_list $NEW_DIR/L1R1_dup.txt

for i in sa pac bwt ann amb; do
cp $OLD_DIR/$FILENAME.$i $NEW_DIR/bwa_index_L1R1/L1R1_dedup.$i
done

################
for d in /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/* ; do
flowcell_id=$(basename $d)

for file in /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/$flowcell_id/*.*; do
name=$(basename $file)
mv $file /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/$flowcell_id/$flowcell_id.$name
done

for file in /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/$flowcell_id/bwa_index_L1R1/*; do
name=$(basename $file)
mv $file /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/$flowcell_id/bwa_index_L1R1/$flowcell_id.$name
done

done




