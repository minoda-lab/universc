args = commandArgs(trailingOnly=TRUE)
print(args)
sam_file <- data.table::fread(as.character(args[1]), data.table = F, skip = 3, header = F, fill=TRUE)
vcf_file <- data.table::fread(as.character(args[2]), data.table = F, skip = 28, fill = TRUE)

head(sam_file)
head(vcf_file)

sam_temp <- sam_file[,12:19]

if(nrow(vcf_file) > 0 ){
# reverse order by position to avoid upstream indesl affecting target position
for(ii in c(nrow(vcf_file)):1){
   print(paste("replace", vcf_file[ii,5], "with", vcf_file[ii,4], "at chr", vcf_file[ii,1], "pos", vcf_file[ii,2]))
   vcf_pos <- vcf_file[ii,2]
   vcf_chr <- vcf_file[ii,1]
   vcf_ref <- vcf_file[ii,4] 
   vcf_alt <- vcf_file[ii,5]
   ref_len <- nchar(vcf_file[ii,4])
   alt_len <- nchar(vcf_file[ii,5])
   sam_pos <- sam_file[,4]
   sam_chr <- sam_file[,3]
   sam_seq <- sam_file[,10]
   sam_qual <- sam_file[,11]
   #print(head(sam_file[,10]))
   read_length <- sapply(sam_file[,10], nchar)
   #print(head(read_length))
   #print(head((sam_chr ==  vcf_chr)))
   #loci <- apply(sam_file, 1, function(read) c(as.numeric(read[4]):as.numeric(read[4])+nchar(read[10])))
#   print(head( sam_pos ))
#   print(head( sam_pos+read_length ))
#   print(head( vcf_pos ))
#   print(head(  (vcf_pos <= sam_pos) & (vcf_pos >= sam_pos+read_length)  ))
#   print(table( (vcf_pos >= sam_pos) & (vcf_pos <= sam_pos+read_length)  ))
   matched_reads <- ((sam_chr ==  vcf_chr) & (vcf_pos >= sam_pos) & (vcf_pos <= sam_pos+read_length))

#   print(which(matched_reads))
#   print(table(matched_reads))
#   print(sam_seq[matched_reads])
#   print(sam_qual[matched_reads])



   print(vcf_pos)
   print(alt_len)
   print(ref_len)
 #  print(sam_pos[matched_reads])
 #  print(sam_pos[matched_reads] + read_length[matched_reads])
 #  print(vcf_pos - sam_pos[matched_reads] + 1)
 #  print( substr(sam_seq[matched_reads], vcf_pos - sam_pos[matched_reads] + 1, vcf_pos - sam_pos[matched_reads] + alt_len) )
 #   print(sam_seq[matched_reads][vcf_pos - sam_pos[matched_reads] + 1])
    print(     substr(sam_seq[matched_reads], vcf_pos - sam_pos[matched_reads] + 1 -1, vcf_pos - sam_pos[matched_reads] + alt_len + 1)   )
 print ( vcf_ref )
   #substr(sam_seq[matched_reads], vcf_pos - sam_pos[matched_reads] + 1, vcf_pos - sam_pos[matched_reads] + alt_len)  <- vcf_ref
   for(jj in c(1:length(sam_seq))[matched_reads]){
         sam_seq[jj] <- paste(c(substr(sam_seq[jj], 0, vcf_pos - sam_pos[jj]),
                        vcf_ref,
                        substr(sam_seq[jj], vcf_pos - sam_pos[jj] + alt_len + 1, nchar(sam_seq[jj]))),
                        collapse = "")
          sam_qual[jj] <- paste(c(substr(sam_qual[jj], 0, vcf_pos - sam_pos[jj]),
                        vcf_ref,
                        substr(sam_qual[jj], vcf_pos - sam_pos[jj] + alt_len + 1, nchar(sam_qual[jj]))),
                        collapse = "")
                 }
     print(     substr(sam_seq[matched_reads], vcf_pos - sam_pos[matched_reads] + 1 -1, vcf_pos - sam_pos[matched_reads] + ref_len + 1)   )
#   print(sam_seq[matched_reads][vcf_pos - sam_pos[matched_reads] + 1])
   #substr(sam_qual[matched_reads], vcf_pos - sam_pos[matched_reads] + 1, vcf_pos - sam_pos[matched_reads] + alt_len) <- paste0(rep("I", ref_len), collapse = "")
 #  print( substr(sam_seq[matched_reads], vcf_pos - sam_pos[matched_reads] + ref_len, vcf_pos - sam_pos[matched_reads] + ref_len ) )
 sam_file[,10] <- sam_seq
 sam_file[,11] <- sam_qual
}

#update CIGAR
read_length <- sapply(sam_file[,10], nchar)
sam_file[,6] <- paste0(read_length, "M")

sam_file[,10] <- sam_seq
sam_file[,11] <- sam_qual

sam_file[sam_file == ""] <- NA
}
#print(sam_file[23,])
#print(sam_file[24,])
sam_file <- sam_file[,1:19]
sam_file[,12:19] <- sam_temp

out_name <- paste0(c(as.character(args[1]), ".converted.sam"), collapse = "")
out_name <- gsub(".sam.converted" ,".converted", out_name)
print(out_name)
write.table(sam_file, file = out_name, col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
