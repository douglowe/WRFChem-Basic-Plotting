;
; Vertical slice plotting script
; 
; This script is for plotting gas-phase chemical data
;
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./vertical_slice_plotting_library.ncl"

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; script settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
root_directory = "/Volumes/Odin/YRD_China/"

input_root_directory  = root_directory+"scenario_output_files/archer_data/base_cases/"
output_root_directory = root_directory+"data_plotting/plots/example_plots/chemistry_vertical_slices/"

scenario = "summer"


plot_type = "pdf"  ; Use this option to output image files
;plot_type = "x11"   ; use this to see the graph (no file will be saved!!!)



;plot_variables = (/ "so4_pm1", "o3", "no", "no2", "hno3", "so2", "ho" /)
;plot_variables = (/ "so4_pm1", "no3_pm1", "nh4_pm1", "nh3", "OM_pm1" /)
;plot_variables = (/ "CLDFRA", "CLDFRA2" /)
;plot_variables = (/ "PM2_5_DRY", "PM10" /)
plot_variables = (/"no2"/)


; This sets the pivot point for the N-S and E-W vertical slices through the domain 
;   The name should be that of the location of the pivot.
;pivot_point = "Shanghai"
;pivot_point@Lat = 31.18
;pivot_point@Lon = 121.45
pivot_point = "Nanjing"
pivot_point@Lat = 32.0
pivot_point@Lon = 118.75

; This sets the height of the vertical cross section (in km)
zmax = 3.0

; Set this to True if you want to use a fixed colour scale for your plotting.
;    Note: If no fixed colour scale is defined for the plotted variable 
;            then a variable colour scale will be used.
fixed_colour_scale = True


; This flag determines if we want to plot a set of defined time points, or the whole campaign.
;    True -> plot the chosen times
;    False -> plot *all* the model data (Warning: it will take 30-50 minutes to load the data)
;    
;    Note: If you set day = "*" then you will plot the selected hour on each day.
;          If you set hour = "*" then you will plot all hours in the selected day. 
;    Note: ***These times are all UTC*** 
plot_selected_times = True
day = (/"13"/)
hour = (/"12"/)
month = (/"*"/)   ; don't change this - always use "*"
year  = (/"*"/)   ; don't change this - always use "*"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; variable settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; set the controls for colour scales here (if required), e.g.:
;;  colour_levels@no2       = (/0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0/)
colour_levels@no2       = (/0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0/)

;;; set the conversion factors for variables here (if required), e.g.:
;; variable_units@T2  = (/"K","C"/)
;;    First string is the original unit (e.g. K - Kelvin)
;;    Second string is the new unit (e.g. C - Celsius)
variable_units@no2 = (/"ppmv","ppbv"/)



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
	plot_resources = True
	load_var_settings(plot_variables(ivar),scenario,\
					var_resources,fixed_colour_scale)

	;;; load variable to plot (converting to the required units)
	plot_var = load_variable(DATA_FILES,plot_variables(ivar),\
								VBS_SPECIES@$plot_variables(ivar)$,variable_units)

	;;; load the grid point height information
	zheight = wrf_user_getvar(DATA_FILES, "z", -1)

	;;; extract vertical slices, and plot the data
	plot_2d_vertical_slices_from_4d_data(DATA_FILES,plot_var,zheight,\
				plot_variables(ivar),output_root_directory+scenario,scenario,\
				var_resources,plot_resources,plot_type,pivot_point,zmax)

end do

end