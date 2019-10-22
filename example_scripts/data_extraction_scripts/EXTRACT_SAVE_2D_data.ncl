;
; 2D data extraction script
; 
; This script is for extracting/generating 2D model data, and saving it to a
;   precreated data file
;
; This script is set to receive optional command line inputs (with default values
;     given below if these are not provided)
;
; The optional inputs are:
;   year, month, day, hour --> strings defining data arrays containing time information
;     (format: "year = (/1997,1998/)" or "year = (/*/)" )
;   input_root_directory, output_root_directory --> strings which define these paths
;   scenario --> string defining the scenario name to use for output subdirectory
;   outfile_name --> string to use for name of the output file
;   domain --> string defining the model domain that is being processed
;
;
; The variables to be created are hard-coded in this script.
;
;


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./data_extraction_library.ncl"

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; script settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
root_directory = "/nerc/n02/n02/lowe/"

if (.not. isvar("input_root_directory")) then
	input_root_directory  = root_directory+"PROMOTE_data_files/Third_VBS_SOA_new_namelist_control/"
end if
if (.not. isvar("output_root_directory")) then
	output_root_directory = root_directory+"PROMOTE_data_analysis/working_space/extracted_data/"
end if

if (.not. isvar("scenario")) then
	scenario = "May_2015-MOZART_boundary_no_dust-base_SOA-BB_Anthro_VBS-x3_emissions-Aging_Rate_1e-11-OC_ratio_tracking"
end if

if (.not. isvar("outfile_name")) then
	outfile_name="test_file.nc"
end if

if (.not. isvar("domain")) then
	domain="d03"
end if


file_variables := (/ "AOD550", "PM1_mass", "OC_ratio", "PM25_mass", "PM10_mass", "PM1_VBS_OM", "PM1_SOA_OM", "PM1_Inorg" /)
proc_variables = (/"AOD550", "DRY_pm1", "ocratio_pm1", "DRY_pm25", "DRY_tot", "POA_pm1", "soa_pm1", "AMSINORG_pm1"/)


; This flag determines if we want to plot a set of defined time points, or the whole campaign.
;    True -> plot the chosen times
;    False -> plot *all* the model data (Warning: it will take 30-50 minutes to load the data)
;    
;    Note: If you set day = "*" then you will plot the selected hour on each day.
;          If you set hour = "*" then you will plot all hours in the selected day. 
;    Note: ***These times are all UTC*** 
plot_selected_times = True
if (.not. isvar("day")) then
	day = (/"*"/)
end if
if (.not. isvar("hour")) then
	hour = (/"00"/)
end if
if (.not. isvar("month")) then
	month = (/"*"/)   ; don't change this - always use "*"
end if
if (.not. isvar("year")) then
	year  = (/"*"/)   ; don't change this - always use "*"
end if

; NOTE: When plotting "rain_hour" you should plot all model data, as this uses the differences
;       between each model output to derive the hourly rain fall.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; variable settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; set the conversion factors for variables here (if required), e.g.:
;; variable_units@T2  = (/"K","C"/)
;;    First string is the original unit (e.g. K - Kelvin)
;;    Second string is the new unit (e.g. C - Celsius)
variable_units@T2  = (/"K","C"/)


;;; set the variable labels here
variable_description@AOD550      = "Total aerosol optical depth in column at 550nm."
variable_unit_description@AOD550 = ""

variable_description@PM1_mass      = "Dry PM1 aerosol mass"
variable_unit_description@PM1_mass = "ug/kg-dryair"

variable_description@OC_ratio      = "PM1 Organic Mass Oxygen:Carbon ratio"
variable_unit_description@OC_ratio = ""

variable_description@PM25_mass   = "Dry PM2.5 aerosol mass"
variable_unit_description@PM25_mass = "ug/kg-dryair"

variable_description@PM10_mass   = "Dry PM10 aerosol mass"
variable_unit_description@PM10_mass = "ug/kg-dryair"

variable_description@PM1_VBS_OM = "PM1 VBS organic mass"
variable_unit_description@PM1_VBS_OM = "ug/kg-dryair"

variable_description@PM1_SOA_OM = "PM1 SOA organic mass"
variable_unit_description@PM1_SOA_OM = "ug/kg-dryair"

variable_description@PM1_Inorg = "PM1 Inorganic mass (without OIN, BC, or Na)"
variable_unit_description@PM1_OM = "ug/kg-dryair"



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; end of settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; check that the output file exists, and load it
output_file = output_root_directory + scenario + "/" + outfile_name
check_file_status(output_file)  
OUT_FILE = addfile(output_file,"w")


;;; open the data files
Dfiles = create_input_file_list(plot_selected_times,domain,year,month,day,hour,input_root_directory+scenario)
check_file_status(Dfiles)  
DATA_FILES = addfiles(Dfiles+".nc","r")

;; copy the Times information to the new file
Times = DATA_FILES[:]->Times
time_dims = dimsizes(Times)
OUT_FILE->Times(0:time_dims(0)-1,:) = (/Times(0:time_dims(0)-1,:)/)


;;; loop through the list of variables that we want to extract
do ivar=0,dimsizes(file_variables)-1
	print("== working on variable "+file_variables(ivar))

	;;; load variable to plot (extracting horizontal slices while we're at it)
	extract_var := load_variable_2D(DATA_FILES,proc_variables(ivar))

	;;; write the data
	if(dimsizes(dimsizes(extract_var)).eq.3)
		OUT_FILE->$file_variables(ivar)$ = (/extract_var/)
	else if(dimsizes(dimsizes(extract_var)).eq.4)
		OUT_FILE->$file_variables(ivar)$ = (/extract_var(:,0,:,:)/)
	end if
	end if
end do

print("SUCCESS: finished processing the data")


end
