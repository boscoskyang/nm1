library(VariantAnnotation)

vcf<-readVcf('C:/Users/syang15/Downloads/GEX1_jeong_normalized.vcf')

if (!is.null(geno(vcf)$GQ)) {
  geno(vcf)$GQ <- matrix(as.numeric(geno(vcf)$GQ) * 100, nrow = nrow(geno(vcf)$GQ))
}

writeVcf(vcf, "C:/Users/syang15/Downloads/GEX1_jeong_normalized_mod.vcf")


#2

vcf<-readVcf('C:/Users/syang15/Downloads/GEX2_han_normalized.vcf')

if (!is.null(geno(vcf)$GQ)) {
  geno(vcf)$GQ <- matrix(as.numeric(geno(vcf)$GQ) * 100, nrow = nrow(geno(vcf)$GQ))
}

writeVcf(vcf, "C:/Users/syang15/Downloads/GEX2_han_normalized_mod.vcf")

#3

vcf<-readVcf('C:/Users/syang15/Downloads/GEX3_gil_normalized.vcf')

if (!is.null(geno(vcf)$GQ)) {
  geno(vcf)$GQ <- matrix(as.numeric(geno(vcf)$GQ) * 100, nrow = nrow(geno(vcf)$GQ))
}

writeVcf(vcf, "C:/Users/syang15/Downloads/GEX3_gil_normalized_mod.vcf")


#4

vcf<-readVcf('C:/Users/syang15/Downloads/GEX4_stroke_normalized.vcf')

if (!is.null(geno(vcf)$GQ)) {
  geno(vcf)$GQ <- matrix(as.numeric(geno(vcf)$GQ) * 100, nrow = nrow(geno(vcf)$GQ))
}

writeVcf(vcf, "C:/Users/syang15/Downloads/GEX4_stroke_normalized_mod.vcf")
