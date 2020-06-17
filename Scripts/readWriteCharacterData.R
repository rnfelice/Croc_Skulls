readCharacterData = function(path, exclude) {
  
  # read the file
  lines = readLines(path)
  
  # find lines corresponding to data
  start_line = grep("MATRIX", toupper(lines), fixed=TRUE) + 1
  end_line   = grep(";", lines, fixed=TRUE)
  end_line   = end_line[which(end_line > start_line)[1]] - 1
  data_lines = lines[start_line:end_line]  
  
  # strip out white space and tabs
  data_lines = gsub("\t"," ", data_lines)
  data_lines = strsplit(data_lines, " ")
  data_lines = lapply(data_lines, function(x) x[x != ""])
  
  # get the taxa
  taxa = sapply(data_lines, function(x) x[1])
  data = do.call(rbind, lapply(data_lines, function(x) as.numeric(x[-1])))
  rownames(data) = taxa

  # excude the specified traits
  data = data[,!1:ncol(data) %in% exclude]
  
  return(data)  
  
}

writeCharacterData <- function(data, file, type) {
  
  lines <- "#NEXUS"
  lines <- c(lines, c(""))
  lines <- c(lines, "BEGIN DATA;")
  lines <- c(lines, paste0("\tDIMENSIONS NTAX=",nrow(data)," NCHAR=",ncol(data),";"))
  if( toupper(type) == "STANDARD") {
    lines <- c(lines, paste0("\tFORMAT DATATYPE=STANDARD SYMBOLS=\"",paste0(sort(unique(as.vector(data))),collapse=""),"\" MISSING=? GAP=-;"))
  } else if ( toupper(type) == "CONTINUOUS" ) {
    lines <- c(lines, paste0("\tFORMAT DATATYPE=CONTINUOUS MISSING=? GAP=-;"))
  } else {
    stop("Invalid data type")
  }
  lines <- c(lines,"MATRIX")
  
  taxa <- rownames(data)
  for(t in taxa) {
    these_traits <- as.numeric(data[t,])
    lines <- c(lines, paste(t, paste0(these_traits, collapse="\t"), sep="\t"))
  }
  
  lines <- c(lines, ";")
  lines <- c(lines, "END;")
  
  # cat(lines, sep="\n")
  
  cat(lines, file=file, sep="\n")
  
}

