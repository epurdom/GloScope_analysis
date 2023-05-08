withCallingHandlers(rmarkdown::render("Code/Process_Datasets/process_perezLupus_default.Rmd",output_file="../../results/Processed_Datasets/process_perezLupus_default.html"),
	error=function(e) print(sys.calls()))
