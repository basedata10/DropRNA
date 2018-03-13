# DropRNA

## Step1: Extract the Cell Barcode

```
python /mnt/gpfs/Database/server/django2/basedata/pipelines/DropSeq/00_get_cellbarcode.py /mnt/gpfs/Database/server/django2/basedata/pipelines/DropSeq/_config.json /mnt/gpfs/Users/zhangxiannian/projects/p11_dropseq/basematics/10X_1.1.fq.gz  ./X10.1.pickle 10X 20 &
```

## Step2: Stats on the Cell Barcode
```
python /mnt/gpfs/Database/server/django2/basedata/pipelines/DropSeq/01_cellbarcode_stats.py _config.json ./X10.1.pickle ./X10.1.stats.json 10X 20000 1000 &
```

## Step3: Split the reads of valid Cell Barcodes
```
nohup python /mnt/gpfs/Database/server/django2/basedata/pipelines/DropSeq/02_split_barcode.py /mnt/gpfs/Users/zhangxiannian/projects/p11_dropseq/basematics/10X_1.1.fq.gz /mnt/gpfs/Users/zhangxiannian/projects/p11_dropseq/basematics/10X_1.2.fq.gz  ./X10.1.stats.json ./reads_10x1 10X &
```

## Step4: Star Alignment

## Step5: Tagging Reads to Gene Names

## Step6: Aggregate the UMI counts and Read counts from a sample to an Matrix
