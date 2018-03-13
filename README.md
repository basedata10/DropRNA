# DropRNA

## Step1: Extract the Cell Barcode

```
python 00_get_cellbarcode.py ./_config.json ./10X_1.1.fq.gz  ./X10.1.pickle 10X 20 &
```

## Step2: Stats on the Cell Barcode
```
python 01_cellbarcode_stats.py _config.json ./X10.1.pickle ./X10.1.stats.json 10X 20000 1000 &
```

## Step3: Split the reads of valid Cell Barcodes
```
python 02_split_barcode.py ./10X_1.1.fq.gz ./10X_1.2.fq.gz  ./X10.1.stats.json ./reads_10x1 10X &
```

## Step4: Star Alignment

## Step5: Tagging Reads to Gene Names

## Step6: Aggregate the UMI counts and Read counts from a sample to an Matrix
