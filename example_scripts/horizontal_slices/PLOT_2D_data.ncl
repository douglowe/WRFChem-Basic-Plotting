;
; 2D data plotting script
; 
; This script is for plotting all single layer data
;
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./horizontal_slice_plotting_library.ncl"

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


;plot_variables = (/ "RAINC", "RAINSH", "RAINNC", "SNOWNC", "GRAUPELNC", "HAILNC", \
;					"rh2", "slp", "td2", "T2", "AOD55", "PBLH" /)

;plot_variables = (/ "RAINC" /)

;plot_variables = (/ "rain_hour" /)

;plot_variables = (/ "OM_tot_belowcloudbase" /)
;plot_variables = (/ "height_belowcloudbase" /)
plot_variables = (/ "rh_belowcloudbase", "depth_belowcloudbase", "height_belowcloudbase", "OM_tot_belowcloudbase" /)



; Set this to True if you want to use a fixed colour scale for your plotting.
;    Note: If no fixed colour scale is defined for the plotted variable 
;            then a variable colour scale will be used.
fixed_colour_scale = True

; Set this to True if you want to plot wind barbs on top of your plot.
wind_switch = True
; This controls the density of the wind barbs (default is 25)
num_barbs = 25


; This flag determines if we want to plot a set of defined time points, or the whole campaign.
;    True -> plot the chosen times
;    False -> plot *all* the model data (Warning: it will take 30-50 minutes to load the data)
;    
;    Note: If you set day = "*" then you will plot the selected hour on each day.
;          If you set hour = "*" then you will plot all hours in the selected day. 
;    Note: ***These times are all UTC*** 
plot_selected_times = True
day = (/"*"/)
hour = (/"00","06","12","18"/)
month = (/"*"/)   ; don't change this - always use "*"
year  = (/"*"/)   ; don't change this - always use "*"

; NOTE: When plotting "rain_hour" you should plot all model data, as this uses the differences
;       between each model output to derive the hourly rain fall.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; variable settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; set the controls for colour scales here (if required), e.g.:
;;  colour_levels@no2       = (/0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0/)
colour_levels@T2 = (/2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40/)
colour_levels@rain_pH  = (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5/)

colour_levels@RHbelowcloudbase = (/0.0,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98/)
colour_levels@rh_belowcloudbase = (/0.0,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98/)

colour_levels@OM_tot_belowcloudbase = (/0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0/)
colour_levels@height_belowcloudbase = (/50,100,200,500,1000,2000,3000,4000,5000,6000,7000/)
colour_levels@depth_belowcloudbase = (/20,50,100,200,400,600,800,1000/)

;;; set the conversion factors for variables here (if required), e.g.:
;; variable_units@T2  = (/"K","C"/)
;;    First string is the original unit (e.g. K - Kelvin)
;;    Second string is the new unit (e.g. C - Celsius)
variable_units@T2  = (/"K","C"/)


;;; set the variable labels here
variable_description@AOD55 = "Total aerosol optical depth in column at 550nm."



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
	map_resources = True
	plot_resources = True
	wind_resources = True
	load_map_settings(plot_variables(ivar),scenario,\
					var_resources,map_resources,wind_resources,\
					num_barbs,fixed_colour_scale)

	;;; load variable to plot (extracting horizontal slices while we're at it)
	plot_var = load_variable_2D(DATA_FILES,plot_variables(ivar))

	;;; plot the data
	plot_2d_data(DATA_FILES,plot_var,\
				plot_variables(ivar),output_root_directory+scenario,scenario,\
				var_resources,map_resources,plot_resources,wind_resources,wind_switch,plot_type)

end do

end
