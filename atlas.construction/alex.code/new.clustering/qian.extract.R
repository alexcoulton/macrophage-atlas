source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(gridExtra)
library(dplyr)
library(Matrix)

g1 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/bc/BC_counts/matrix.mtx'
g2 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/bc/BC_counts/genes.tsv'
g3 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/bc/BC_counts/barcodes.tsv'

qian.bc = ReadMtx(mtx = g1, cells = g3, features = g2)

g1 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/crc/CRC_counts/matrix.mtx'
g2 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/crc/CRC_counts/genes.tsv'
g3 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/crc/CRC_counts/barcodes.tsv'

qian.crc = ReadMtx(mtx = g1, cells = g3, features = g2)

g1 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/ovc/OvC_counts/matrix.mtx'
g2 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/ovc/OvC_counts/genes.tsv'
g3 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/ovc/OvC_counts/barcodes.tsv'

qian.ovc = ReadMtx(mtx = g1, cells = g3, features = g2)

g1 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/lung/LC_counts/matrix.mtx'
g2 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/lung/LC_counts/genes.tsv'
g3 = '~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/qian/lung/LC_counts/barcodes.tsv'

qian.lung = ReadMtx(mtx = g1, cells = g3, features = g2)


dim(qian.bc)
dim(qian.crc)
dim(qian.ovc)
dim(qian.lung)


