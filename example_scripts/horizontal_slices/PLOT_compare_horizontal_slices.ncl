;
; Horizontal slice plotting script
; 
; This script is for plotting gas-phase chemical data
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



;plot_variables = (/ "so4_pm1", "o3", "no", "no2", "hno3", "so2", "ho" /)
;plot_variables = (/ "so4_pm1", "no3_pm1", "nh4_pm1", "nh3", "OM_pm1" /)
;plot_variables = (/ "CLDFRA", "CLDFRA2" /)
;plot_variables = (/ "PM2_5_DRY", "PM10" /)
;plot_variables = (/"oin_tot","so4_tot","OM_tot"/)
plot_variables = (/"no"/)


; This sets the altitude levels at which you wish to extract the model data
;   Note 1: if you want the data at ground level, then use an altitude of 0.0
;   Note 2: these heights are in "metres above sea level"
;altitude_levels = (/ 0.0, 50.0, 100.0, 300.0, 500.0 /)
altitude_levels = (/ 0.0 /)

; Set this to True if you want to use a fixed colour scale for your plotting.
;    Note: If no fixed colour scale is defined for the plotted variable 
;            then a variable colour scale will be used.
fixed_colour_scale = True

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
day = (/"1*"/)
hour = (/"00","06","12","18"/)
month = (/"*"/)   ; don't change this - always use "*"
year  = (/"*"/)   ; don't change this - always use "*"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; variable settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; set the controls for colour scales here (if required), e.g.:
;;  colour_levels@no2       = (/0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0/)
colour_levels@no2       = (/-50.0,-20.0,-10.0,-5.0,-2.0,-1.0,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0/)
colour_levels@o3       = (/-275,-225,-175,-125,-75,-25,25,75,125,175,225,275/)

;;; set the conversion factors for variables here (if required), e.g.:
;; variable_units@T2  = (/"K","C"/)
;;    First string is the original unit (e.g. K - Kelvin)
;;    Second string is the new unit (e.g. C - Celsius)
variable_units@no2 = (/"ppmv","ppbv"/)
variable_units@o3 = (/"ppmv","ppbv"/)


;;; set the variable labels here
variable_description@ho    = "OH"
variable_description@o3    = "O~B~3~N~"
variable_description@no    = "NO"
variable_description@no2   = "NO~B~2~N~"
variable_description@TAU55 = "Aerosol Optical Thickness at 550nm"


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

	;;; load variables to plot (extracting horizontal slices while we're at it)
	plot_var_a = load_variable_and_sample_at_given_altitudes(DATA_FILES_A,plot_variables(ivar),altitude_levels)
	plot_var_b = load_variable_and_sample_at_given_altitudes(DATA_FILES_B,plot_variables(ivar),altitude_levels)

	;;; get the difference between the two data sets (A - B)
	plot_var = difference_calculation(plot_var_a,plot_var_b,percentage_difference)
	

	;;; plot the data
	plot_2d_horizontal_slice_from_4d_data(DATA_FILES_A,plot_var,\
				plot_variables(ivar),output_root_directory+scenario,scenario,\
				var_resources,map_resources,plot_resources,False,0,plot_type,altitude_levels)

end do

end
