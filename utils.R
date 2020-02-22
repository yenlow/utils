# Title     : utilities to partially source R scripts by regex/lines, save graphics etc
# Created by: yensia-low
# Created on: 07-Feb-20

# https://stackoverflow.com/questions/26245554/execute-a-set-of-lines-from-another-r-file
sourcePartial <- function(file,startTag='#from here',endTag=NULL, startskip=0, endskip=0) {
  lines = scan(file, what=character(), sep="\n", quiet=TRUE)
  st = grep(startTag,lines)
  en = ifelse(is.null(endTag), length(lines), grep(endTag,lines))
  tc = textConnection(lines[(st+startskip):(en-endskip)])
  source(tc)
  close(tc)
}

saveGraphics <- function(file, outdir="output/plots",
                         startTag='^png\\(|^pdf\\(',
                         endTag="^dev\\.off\\(\\)$") {
  lines = scan(file, what=character(), sep="\n", quiet=TRUE)

  #find the lines for saving graphics using the graphics device
  st = grep(startTag,lines)
  en = grep(endTag,lines)
  indices = data.frame(st,en)

  graphicsLinesAll = c()
  outpaths = c()
  for(i in 1:nrow(indices)){
    ##gsub the original path with outdir
    replaceStr = paste0('("', outdir,'/')
    lineStart = gsub('\\(\\"[[:alnum:]_]+/', replaceStr, lines[indices$st[i]])
    outpath  = regmatches(lineStart,regexpr('\\"(.*?)\\"',lineStart))
    #explicitly print the object so they will be rendered and saved
    lineb4devoff = paste0("print(",lines[indices$en[i]-1],")")

    #these will the be new lines that will be sourced
    graphicsLines = lines[indices$st[i]:indices$en[i]]
    graphicsLines[1] = lineStart
    graphicsLines[length(graphicsLines)-1] = lineb4devoff

    graphicsLinesAll = c(graphicsLinesAll, graphicsLines)
    outpaths = c(outpaths, outpath)
  }

    tc = textConnection(graphicsLinesAll)
    source(tc)
    close(tc)
    print("Graphics were saved to:")
    print(outpaths)
}

# Returns middle of a range, e.g. 21-30 => 25
mid <- function(x) floor(mean(as.numeric(x)))

ageChecker <- function(v){
  #strip white spaces
  v = gsub("\\s", "", v)

  #print frequency of invalid ages
  badAge = grep("(^\\d{1,3}$)", v, value=TRUE, invert=TRUE)
  print("Invalid ages:")
  print(table(badAge,useNA="ifany"))

  #change ranges to medians
  ageRangeID = which(grepl("\\d+-\\d+", v))
  v[ageRangeID] = as.numeric(lapply(strsplit(v[ageRangeID],'-'),mid))

  #change 20s to 25
  v = gsub("\\ds$", "5", v)
  #drop non-numbers or empty strings
  v = gsub("^\\D+$|''", NA, v)
  v = floor(as.numeric(v))

  #print cleaned up dates
  print("Cleaned up ages:")
  print(table(v,useNA="ifany"))
}

