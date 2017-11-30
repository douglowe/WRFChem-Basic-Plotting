;
;  Plot domain statistical data as box and whisker plots 
;
;  The box middle point is the "median" value.
;  The box upper and lower edges are the 75th and 25th percentile values
;  The upper and lower whiskers are the 95th and 5th percentile values
;


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./time_series_plotting_library.ncl"

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; script settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
root_directory = "/nerc/n02/n02/lowe/"

input_root_directory  = root_directory+"RONOCO_VBS_data_files/"
output_root_directory = root_directory+"RONOCO_VBS_data_analysis/working_space/plots/example_plots/chemistry_timeseries/"

scenario = "wet_data"


plot_type = "pdf"  ; Use this option to output image files
;plot_type = "x11"   ; use this to see the graph (no file will be saved!!!)



;plot_variables = (/ "o3", "no", "no2", "hno3", "so2", "ho" /)
plot_variables = (/ "RHbelowcloudbase" /)
;plot_variables = (/ "so4_pm1", "no3_pm1", "nh4_pm1", "nh3", "OM_pm1" /)
;plot_variables = (/ "CLDFRA", "CLDFRA2" /)
;plot_variables = (/ "PM2_5_DRY", "PM10" /)
;plot_variables = (/"oin_pm25"/)
;plot_variables = (/"TAU55"/)
;plot_variables = (/"no"/)


; This sets the altitude levels at which you wish to extract the model data
;   Note 1: if you want the data at ground level, then use an altitude of 0.0
;   Note 2: these heights are in "metres above sea level"
;altitude_levels = (/ 0.0, 50.0, 100.0, 300.0, 500.0, 1000.0, 2000.0 /)
;altitude_levels = (/ 0.0 /)
altitude_levels = (/ 0.0, 500.0 /)

; Set this to True if you want to "trim" the edge of the domain before calculating 
;    the domain stats.
;  "trim" is the number of grid cells on each side of the domain to remove.
domain_trim = True
domain_trim@trim = 5

; Set this to True if you want to set the y-axis limits below 
;
yaxis_limits = True

; Set this to True if you want to label plots with local time, not UTC
;    You will need to set the offset to do this. 
use_local_time = False
use_local_time@offset = 8




; This flag determines if we want to plot a set of defined time points, or the whole campaign.
;    True -> plot the chosen times
;    False -> plot *all* the model data (Warning: it will take 30-50 minutes to load the data)
;    
;    Note: If you set day = "*" then you will plot the selected hour on each day.
;          If you set hour = "*" then you will plot all hours in the selected day.
;    Note: ***These times are all UTC*** 
plot_selected_times = True
day = (/"*"/)
hour = (/"00","08","16"/)
month = (/"*"/)   ; don't change this - always use "*"
year  = (/"*"/)   ; don't change this - always use "*"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; variable settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; set y-axis limits here (if required), e.g.:
;;  yaxis_limits@no2       = (/0.0,100.0/)
yaxis_limits@o3       = (/-1.0,100.0/)
yaxis_limits@RHbelowcloudbase       = (/80.0,101.0/)


;;; set the conversion factors for variables here (if required), e.g.:
;; variable_units@T2  = (/"K","C"/)
;;    First string is the original unit (e.g. K - Kelvin)
;;    Second string is the new unit (e.g. C - Celsius)
variable_units@no2   = (/"ppmv","ppbv"/)
variable_units@no    = (/"ppmv","ppbv"/)
variable_units@o3    = (/"ppmv","ppbv"/)
variable_units@no3   = (/"ppmv","pptv"/)
variable_units@T2    = (/"K","C"/)


;;; set the variable labels here
variable_description@ho    = "OH"
variable_description@o3    = "O~B~3~N~"
variable_description@no    = "NO"
variable_description@no2   = "NO~B~2~N~"
variable_description@TAU55 = "Aerosol Optical Thickness at 550nm"



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; end of settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; check that the output directory exists, if not then create it
system("if ! test -d  "+output_root_directory+" ; then mkdir "+output_root_directory+"; fi")
system("if ! test -d  "+output_root_directory+scenario+" ; then mkdir "+output_root_directory+scenario+"; fi")



;;; open the data files
Dfiles = create_input_file_list(plot_selected_times,"d01",year,month,day,hour,input_root_directory+scenario)
check_file_status(Dfiles)  
DATA_FILES = addfiles(Dfiles+".nc","r")


;;; loop through the list of variables that we want to plot
do ivar=0,dimsizes(plot_variables)-1
	print("== working on variable "+plot_variables(ivar))

	;;; load the settings for this variable
	var_resources = True
	load_plot_settings(plot_variables(ivar),scenario,var_resources,use_local_time)

	;;; load variable to plot (extracting horizontal slices while we're at it)
	plot_var = load_2D_3D_variable_and_sample_at_given_altitudes(DATA_FILES,plot_variables(ivar),altitude_levels)

	;;; calculate the domain stats
	plot_stats = calculate_domain_stats(plot_var,domain_trim)

	;;; plot the data
	plot_stats_data(DATA_FILES,plot_stats,\
				plot_variables(ivar),output_root_directory+scenario,scenario,\
				var_resources,plot_type,altitude_levels,yaxis_limits,use_local_time)
				
	;;; write out the plotted data as a text file
	write_stats_data(DATA_FILES,plot_stats,\
				plot_variables(ivar),output_root_directory+scenario,scenario,\
				altitude_levels,use_local_time)

	;;; tidy up
	delete([/plot_var,plot_stats/])

end do

end
