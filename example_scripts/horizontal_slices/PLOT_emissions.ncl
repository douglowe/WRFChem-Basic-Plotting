;
; Script for plotting the emission data stored in the wrfout files.
; 
; Use this script to double-check the emissions that you have used in
; each scenario are what you intended to use.
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./horizontal_slice_plotting_library.ncl"

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; script settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
root_directory = "/Volumes/Odin/YRD_China/"

;input_root_directory  = root_directory+"scenario_output_files/base_cases/"
input_root_directory = root_directory+"scenario_output_files/archer_data/zhenze/"
output_root_directory = root_directory+"data_plotting/plots/example_plots/chemistry_horizontal_slices/"

;scenario = "winter"
scenario = "simulation2"


;plot_type = "pdf"  ; Use this option to output image files
plot_type = "x11"   ; use this to see the graph (no file will be saved!!!)



plot_variables = (/"E_OC"/)


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
day = (/"14"/)
hour = (/"00"/)
month = (/"*"/)   ; don't change this - always use "*"
year  = (/"*"/)   ; don't change this - always use "*"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; variable settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; set the controls for colour scales here (if required), e.g.:
;;  colour_levels@no2       = (/0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0/)
colour_levels@E_OC       = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5/)
colour_levels@E_OC_DOM   = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5/)
colour_levels@E_OC_TRA   = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5/)
colour_levels@E_OC_25_10 = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
colour_levels@E_ALD      = (/0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6/)
colour_levels@E_BC_1     = (/0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1/)
colour_levels@E_C2H5OH   = (/2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0/)
colour_levels@E_CH3OH    = (/2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0/)
colour_levels@E_CO       = (/500.0,1000.0,1500.0,2000.0,2500.0,3000.0,3500.0,4000.0,4500.0,5000.0,5500.0,6000.0,6500.0,7000.0,7500.0,8000.0/)
colour_levels@E_CSL      = (/0.0,1.0/)
colour_levels@E_EC_1_25  = (/0.0,1.0/)
colour_levels@E_EC_25_10 = (/0.0,1.0/)
colour_levels@E_ETH      = (/2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0/)
colour_levels@E_HC3      = (/0.0,1.0/)
colour_levels@E_HC5      = (/0.0,1.0/)
colour_levels@E_HC8      = (/0.0,1.0/)
colour_levels@E_HCHO     = (/0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0,3.3,3.6,3.9,4.2,4.5,4.8,5.1,5.4,5.7,6.0/)
colour_levels@E_ISO      = (/0.0,1.0/)
colour_levels@E_KET      = (/5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0/)
colour_levels@E_NH3      = (/5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0/)
colour_levels@E_NO2      = (/10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0/)
colour_levels@E_NO       = (/100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0/)
colour_levels@E_OIN_10   = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5/)
colour_levels@E_OIN_25   = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5/)
colour_levels@E_OL2      = (/5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0/)
colour_levels@E_OLI      = (/2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0/)
colour_levels@E_OLT      = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0/)
colour_levels@E_ORA2     = (/0.0,1.0/)
colour_levels@E_PAR      = (/100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1100.0,1200.0/)
colour_levels@E_PM25     = (/0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6/)
colour_levels@E_SO2      = (/50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0,550.0,600.0/)
colour_levels@E_TOL      = (/0.0,1.0/)
colour_levels@E_XYL      = (/10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0/)

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
	map_resources = True
	plot_resources = True
	wind_resources = True
	load_map_settings(plot_variables(ivar),scenario,\
					var_resources,map_resources,wind_resources,\
					0,fixed_colour_scale)

	;;; load variable to plot (extracting horizontal slices while we're at it)
	plot_var = load_emissions(DATA_FILES,plot_variables(ivar))
	
	;;; plot the data
	plot_2d_horizontal_slice_from_4d_data(DATA_FILES,plot_var,\
				plot_variables(ivar),output_root_directory+scenario,scenario,\
				var_resources,map_resources,plot_resources,wind_resources,False,\
				plot_type,0.0)

end do

end