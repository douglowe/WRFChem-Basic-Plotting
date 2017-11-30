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
root_directory = "/Volumes/Odin/YRD_China/"

input_root_directory_a = root_directory+"scenario_output_files/archer_data/base_cases/"
input_root_directory_b = root_directory+"scenario_output_files/archer_data/zhenze/"
output_root_directory  = root_directory+"data_plotting/plots/example_plots/chemistry_horizontal_slices/"

scenario_a = "summer"
scenario_b = "simulation3"
scenario = scenario_b+"-"+scenario_a


;plot_type = "pdf"  ; Use this option to output image files
plot_type = "x11"   ; use this to see the graph (no file will be saved!!!)


;plot_variables = (/ "RAINC", "RAINSH", "RAINNC", "SNOWNC", "GRAUPELNC", "HAILNC", \
;					"rh2", "slp", "td2", "T2", "AOD55", "PBLH" /)

;plot_variables = (/ "RAINC" /)

;plot_variables = (/ "rain_hour" /)

plot_variables = (/ "T2","PBLH","rain_hour" /)

; Set this to True if you want to use a fixed colour scale for your plotting.
;    Note: If no fixed colour scale is defined for the plotted variable 
;            then a variable colour scale will be used.
fixed_colour_scale = False

; Set this to True if you want to plot the percentage difference, instead of the absolute difference.
;
percentage_difference = True


; This flag determines if we want to plot a set of defined time points, or the whole campaign.
;    True -> plot the chosen times
;    False -> plot *all* the model data (Warning: it will take 30-50 minutes to load the data)
;    
;    Note: If you set day = "*" then you will plot the selected hour on each day.
;          If you set hour = "*" then you will plot all hours in the selected day. 
;    Note: ***These times are all UTC*** 
plot_selected_times = True
day = (/"13"/)
hour = (/"11","12"/)
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

;;; this color map is suitable for difference plots!
colormap = "testcmap"




;;; check that the output directory exists, if not then create it
system("if ! test -d  "+output_root_directory+" ; then mkdir "+output_root_directory+"; fi")
system("if ! test -d  "+output_root_directory+scenario+" ; then mkdir "+output_root_directory+scenario+"; fi")



;;; open the data files
Dfiles_a = create_input_file_list(plot_selected_times,"d01",year,month,day,hour,input_root_directory_a+scenario_a)
Dfiles_b = create_input_file_list(plot_selected_times,"d01",year,month,day,hour,input_root_directory_b+scenario_b)
check_file_status(Dfiles_a)  
check_file_status(Dfiles_b)  
DATA_FILES_A = addfiles(Dfiles_a+".nc","r")
DATA_FILES_B = addfiles(Dfiles_b+".nc","r")


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
					0,fixed_colour_scale)

	;;; load variable to plot (extracting horizontal slices while we're at it)
	plot_var_a = load_variable_2D(DATA_FILES_A,plot_variables(ivar))
	plot_var_b = load_variable_2D(DATA_FILES_B,plot_variables(ivar))

	;;; get the difference between the two data sets (A - B)
	plot_var = difference_calculation(plot_var_a,plot_var_b,percentage_difference)


	;;; plot the data
	plot_2d_data(DATA_FILES_A,plot_var,\
				plot_variables(ivar),output_root_directory+scenario,scenario,\
				var_resources,map_resources,plot_resources,False,0,plot_type)

end do

end
