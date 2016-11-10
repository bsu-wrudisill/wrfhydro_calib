## JLM: dumping code to be added/edited for 2D processes at the top
## gather all obs before hand

#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

.libPaths("/glade/u/home/adugger/system/R/Libraries/R3.2.2")
library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)

#########################################################
# SETUP
#########################################################

source("calib_utils.R")

### WILL: sets up parallel processing

# Multi-core
parallelFlag <- TRUE
ncores <- 8
if (parallelFlag && ncores>1) {
        library(doParallel)
        cl <- makeForkCluster(ncores)
        registerDoParallel(cl)
}


#########################################################
                            # MAIN CODE
#########################################################

### WILL: proj_data.Rdata == where previous calib outputs are stored?

# First loop check
if (file.exists("proj_data.Rdata")) { 
   load("proj_data.Rdata")
} else {
   # First run so need to initialize
   ReadNamelist("namelist.calib")
   cyclecount <- 0

   # Setup plot directory
   writePlotDir <- paste0(runDir, "/plots")
   dir.create(writePlotDir)

   ## WILL: What goes into obs file? I'm guessing 'observations'... Yes. 

   # Load obs so we have them for next iteration
   load(obsFile)
   #obsDT$q_cms <- NULL


###### WILL: This part will be removed/changed
######
   # Find the index of the gage
  # rtLink <- ReadRouteLink(rtlinkFile)
  # rtLink <- data.table(rtLink)
  # linkId <- which(trimws(rtLink$gages) %in% siteId)
##########
#########
   
   ####WILL: this chunk should stay the same. 'paramBnds' must come from
   #### the calib.utils.R script... 

   # Setup value lists from paramBnds
   xnames <- paramBnds$param
   x0 <- paramBnds$ini
   names(x0) <- xnames
   x_min <- paramBnds$min
   names(x_min) <- xnames
   x_max <- paramBnds$max
   names(x_max) <- xnames

   #WILL: This dataframe contains all of the iterations of Parameters

   # Initialize parameter archive DF
   message("Initialize parameter archive")
   x_archive <- as.data.frame(matrix(, nrow=1, ncol=length(xnames)+2))
   names(x_archive) <- c("id", xnames, "obj")

   # Output parameter set
   x_new <- x0
   cyclecount <- 1

   x_new_out <- c(cyclecount, x_new)
   names(x_new_out)[1] <- "id"
   write.table(data.frame(t(x_new_out)), file="params_new.txt", row.names=FALSE, sep=" ")

   dir.create("archive")

   # Save and exit
   save.image("proj_data.Rdata")
   quit("no")
}


if (cyclecount > 0) {

   file.rename("params_new.txt", paste0("CALIB_RESULTS/OUTPUT", cyclecount, "/params_new.txt"))

   # Read model out and calculate performance metric
   outPath <- paste0(runDir, "/CALIB_RESULTS/OUTPUT", cyclecount)
   print(outPath)

   # Read files
   message("Reading model out files.")

########
########

   system.time({
   
   #collect observation files and reduce name times to the LDASOUT convention

   obsDir <- "~/tlmod"
   oFiles <- list.files(obsDir, patt='MYD10A1*')

   #for the synthetic data case, obsTimePatterns <- names(oList)
   
   #specific case for modis data.. maybe the best thing to to is rename to Obs files beforehand 
   #to the same specification as model out files, i.e yyyymmddhhmm.filename
   oList <- paste0(as.Date(as.integer(substr(basename(oFiles),14,16)), origin=paste0(substr(basename(oFiles), 10,13), "-01-01")),"_12:00")
   obsTimePatterns <- gsub('(-|_| |:)', '' , oList)

   ##this way time and the file name is matched. maybe change var name...
   obsTimePatterns <- paste0("(",paste0(paste0(obsTimePatterns,'.LDASOUT_DOMAIN1', collaps=''),collapse='|'),")")

   modelFiles <- list.files(path = outputPath, 
                            patt=obsTimePatterns, 
                            full=TRUE)

   if (length(modelFiles) == 0) stop("No matching files in specified directory.")
   
   allLdasOut <- plyr::llply(modelFiles, ReadLdasOut)
   
   ####not very elegant
   names(allLdasOut)<-oList



########
########

   #####WILL: Should spit out dates from the outfile names...
   #filesListDate <- as.POSIXct(unlist(plyr::llply(strsplit(basename(filesList),"[.]"), '[',1)), format = "%Y%m%d%H%M", tz = "UTC")
   
   ####WIll: this allows for model spinup, outfiles aren't read from before startdate(?)
   #whFiles <- which(filesListDate >= startDate)

   ###WILL Subsets outfiles that are newer than start date
   #filesList <- filesList[whFiles]
   #if (length(filesList) == 0) stop("No matching files in specified directory.")


#######WILL: this will need to change... ReadChFile Refers to the util function which reads 
#######streamflow from the output Netcdf file.  

   #chrt <- as.data.table(plyr::ldply(filesList, ReadChFile, linkId, .parallel = parallelFlag))
   #})

   #######WILL: This will also need to change
   # Convert to daily
   #chrt.d <- Convert2Daily(chrt)
   #chrt.d[, site_no := siteId]

   ###WILL: not sure what assign is doing
   assign(paste0("chrt.d.", cyclecount), chrt.d)
   
   ###WILL: saves the chrt (the output discharge) to .Rdata file
   save(list=c(paste0("chrt.d.", cyclecount)), file=paste0("archive/", paste0("chrt.d.", cyclecount), ".Rdata"))

   # Merge
   ###WILL: where does setkey come from? 
   setkey(chrt.d, "site_no", "POSIXct")
   setkey(obsDT, "site_no", "POSIXct")

   ##WILL: 'merge' merges two dataframes by parameter. This step puts obs and output 
   ## the same dataframe..
   chrt.d <- merge(chrt.d, obsDT, by=c("site_no", "POSIXct"), all.x=FALSE, all.y=FALSE)
   
   ##WILL: objFn (note: not a verb) returns the objective function
   F_new <- objFn(chrt.d$q_cms, chrt.d$obs)
   print(F_new)



#############################################################################
#######         WILL: Below should not need to change       #################

   # Archive results
   x_archive[cyclecount,] <- c(cyclecount, x_new, F_new)

   # Evaluate performance metric
   # WILL: recall, x_new is the vector of params
   if (cyclecount == 1) {
      x_best <- x_new
      F_best <- F_new
   } else if (F_new <= F_best) {
      x_best <- x_new
      F_best <- F_new
   }

   ##WILL: where does m come from?
   if (cyclecount < m) {
      # Select next parameter set
      
      ###WILL: DDS is the algorithm that takes the best, min, and max and finds params for 
      ###the next timestep. where does r come from?

      x_new <- DDS.sel(i=cyclecount, m=m, r=r, xnames=xnames, x_min=x_min, x_max=x_max, x_best=x_best)
      cyclecount <- cyclecount+1  

      # Output next parameter set
      x_new_out <- c(cyclecount, x_new)
      names(x_new_out)[1] <- "id"
      write.table(data.frame(t(x_new_out)), file="params_new.txt", row.names=FALSE, sep=" ")
   }

   # Stop cluster
   if (parallelFlag) stopCluster(cl)

   # Update plot
   gg <- ggplot(data=x_archive, aes(x=id, y=obj)) + 
              geom_point() + theme_bw() + 
              labs(x="run", y="objective function")
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_calib_run_obj.png"),
              plot=gg, units="in", width=6, height=5, dpi=300)

   # Archive output
   if (!archiveOutput) system(paste0("rm -r ", outPath), intern=FALSE)

   # Archive model run dir
   modFromPath <- paste0(runDir, "/RUN.CALTMP")
   modToPath <- paste0(runDir, "/CALIB_RUNS/RUN.CALTMP", cyclecount)
   if (archiveRun) {
      system(paste0("mv ", modFromPath, " ", modToPath), intern=FALSE)
   } else {
      system(paste0("rm -r ", modFromPath), intern=FALSE)
   }

   # Save and exit
   ###WILL: what gets saved in 'save.image'?....

   save.image("proj_data.Rdata")
   quit("no")

}



