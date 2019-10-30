;
; library file for the scripts which deal with horizontal plotting of model data.
;
;  
;
;
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

load "../../library/variable_loading_library.ncl"

begin

;;;;;;; Setting up the logical variables which we will use to carry 
;;;;;;;   attributes controlling the variable setups
;; used for carrying colour scale information, e.g.:
;;   colour_levels@no2       = (/0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0/) 
colour_levels = True

;; used for carrying information about changing variable units, e.g.:
;;     variable_units@T2  = (/"K","C"/)
variable_units = True

;; used for carrying the variable descriptions
variable_description = True

;; used for carrying the variable unit description
variable_unit_description = True


;;; define if the species is a VBS compound or not
VBS_SPECIES = True
VBS_SPECIES@OM_pm1  = True
VBS_SPECIES@OM_pm25 = True
VBS_SPECIES@OM_tot  = True
VBS_SPECIES@OM_tot_belowcloudbase = True
VBS_SPECIES@ALL_pm1 = True
VBS_SPECIES@DRY_pm1 = True
VBS_SPECIES@ocratio_pm1 = True
VBS_SPECIES@ioratio_pm1 = True
VBS_SPECIES@amsratio_pm1 = True
VBS_SPECIES@ALL_pm25 = True
VBS_SPECIES@DRY_pm25 = True
VBS_SPECIES@DRY_tot = True


;;; default colormap - this could be changed in the calling script if needed
colormap = "WhBlGrYeRe"




;;;;;;; settings loading routine
undef("load_map_settings")
procedure load_map_settings(var:string,scenario:string,\
					resources:logical,map_resources:logical,wind_resources:logical,\
					num_barbs:numeric,fixed_colour_scale:logical) 
begin
	
	;;; delete the previous resources
	delete_VarAtts(resources,-1)
	
	;;; set general plotting resources
	resources@MainTitle 	= ""
	resources@InitTime 		= False
	resources@Footer 		= False
	resources@cnFillOn 		= True
	resources@cnFillMode	= "RasterFill"
	
	if(fixed_colour_scale .and. isatt(colour_levels,var))
		resources@cnLevelSelectionMode = "ExplicitLevels"
		resources@cnLevels = colour_levels@$var$
	end if

	;;; map settings
	map_resources@mpGeophysicalLineColor = "Black"
	map_resources@mpNationalLineColor    = "Black"
	map_resources@mpGeophysicalLineThicknessF = 8.0
	
	;;; wind settings
	wind_resources@NumVectors = num_barbs

end







;;;;;;; variable loading routines
undef("load_emissions")
function load_emissions(FILES:list,var:string)
local var1, var2, var3, variable
begin


	if (var .eq. "E_OC")
		;;; load data
		var1 = FILES[:]->E_OC_DOM
		var2 = FILES[:]->E_OC_TRA
		var3 = FILES[:]->E_OC_25_10
		variable = var1 + var2 + var3
		copy_VarAtts(var1,variable)
	
		variable@description = "EMISSION OF ORGANIC CARBON"

	else ;; otherwise just load generic data

		;;; load data
		variable = FILES[:]->$var$

	end if
            

	return(variable)	
	
end





undef("load_variable_and_sample_at_given_altitudes")
function load_variable_and_sample_at_given_altitudes(FILES:list,var:string,altitude_levels:numeric)
local var1, var2, var3, var_temp, variable, alt_km
begin

	;;; check if we are processing a VBS variable
	if(isatt(VBS_SPECIES,var))
		;;; call routine for extracting whole domain volume
		var_temp = load_variable(FILES,var,VBS_SPECIES@$var$,variable_units)		
	else
		;;; call routine for extracting whole domain volume, assume it's not a VBS species
		var_temp = load_variable(FILES,var,False,variable_units)
	end if

	;;; convert from metres to km for the altitude
	alt_km = altitude_levels / 1000.0


	;;; extract horizontal slices from the domain volume
	variable = wrf_user_vert_interp(FILES,var_temp,"ght_msl",alt_km,True)

	;;; replace any altitudes of zero with the data from the ground level
	do in = 0,dimsizes(alt_km)-1
		if(alt_km(in) .eq. 0.0)
			variable(:,in,:,:) = (/var_temp(:,0,:,:)/)
		end if
	end do



	return(variable)	
	
end



undef("load_variable_2D")
function load_variable_2D(FILES:list,var:string)
local var1, var2, var3, var4, var5, var6, vara, dima, itime, var_temp, vbs_flag
begin


	;; obtain a rough rain pH value
	if(var .eq. "rain_pH")

		var_temp = get_rain_approx_ph(FILES,-1)
	


	;; process the 3D aerosol optical thickness data to give a total column aerosol optical depth
	else if(var .eq. "AOD55")
		var1 = load_variable(FILES,"TAU55",False,variable_units)
	
		; Find total column depth by summing through vertical levels:
		; declare variable to hold column data for variable. Set to ground level first
		var_temp = var1(:,0,:,:)	
		
		var_temp@_FillValue = -9999.
					
		dima = dimsizes(var1)
		
		do ik = 1, dima(1) - 1	; loop through vertical levels
			var_temp = var_temp + var1(:,ik,:,:)
		end do
		

	
	;; process the accumulated precipitation data to give hourly precip rates
	else if(isStrSubset(var,"_hour"))
		if(var .eq. "RAINC_hour")
			vara = load_variable(FILES,"RAINC",False,variable_units)
			vara@description = "TOTAL CUMULUS PRECIPITATION"	
		else if(var .eq. "RAINSH_hour")
			vara = load_variable(FILES,"RAINSH",False,variable_units)
			vara@description = "SHALLOW CUMULUS PRECIPITATION"	
		else if(var .eq. "RAINNC_hour")
			vara = load_variable(FILES,"RAINNC",False,variable_units)
			vara@description = "TOTAL GRID SCALE PRECIPITATION"	
		else if(var .eq. "SNOWNC_hour")
			vara = load_variable(FILES,"SNOWNC",False,variable_units)
			vara@description = "TOTAL GRID SCALE SNOW AND ICE"	
		else if(var .eq. "GRAUPELNC_hour")
			vara = load_variable(FILES,"GRAUPELNC",False,variable_units)
			vara@description = "TOTAL GRID SCALE GRAUPEL"	
		else if(var .eq. "HAILNC_hour")
			vara = load_variable(FILES,"HAILNC",False,variable_units)
			vara@description = "TOTAL GRID SCALE HAIL"	
		else if(var .eq. "rain_hour" .or. var .eq. "precip_hour")
			var1 = load_variable(FILES,"RAINC",False,variable_units)	
			var2 = load_variable(FILES,"RAINSH",False,variable_units)	
			var3 = load_variable(FILES,"RAINNC",False,variable_units)	
			
			if(var .eq. "precip_hour")
				var4 = load_variable(FILES,"SNOWNC",False,variable_units)	
				var5 = load_variable(FILES,"GRAUPELNC",False,variable_units)	
				var6 = load_variable(FILES,"HAILNC",False,variable_units)	
				vara = var1+var2+var3+var4+var5+var6
				vara@description = "Total Precipitation"
			else	
				vara = var1+var2+var3
				vara_description = "Total Rain Precipitation"
			end if
		else
			print("species "+var+" is unrecognised")
			break
		end if
		end if
		end if
		end if
		end if
		end if
		end if
	
		;; get the dimensions of the data set
		dima = dimsizes(vara)
	
		;; initialise the output variable
		var_temp = vara
		;; loop through the output variable, calculating the hourly precip rate
		do itime= 1,dima(0)-1
			var_temp(itime,:,:) = vara(itime,:,:) - vara(itime-1,:,:)
		end do
		;; set the first value of the output variable to zero
		var_temp(0,:,:) = 0.0 
		
		var_temp@units = "mm/hr"
	
		
	else ;;; extracting all other 2D datasets
	
		
		if(isatt(VBS_SPECIES,var))
			vbs_flag = VBS_SPECIES@$var$
		else
			vbs_flag = False
		end if
	
		;;; call routine for extracting whole domain volume, assume it's not a VBS species
		var_temp = load_variable(FILES,var,vbs_flag,variable_units)
	
	end if
	end if
	end if


	return(var_temp)


end


;;;;;;; data writing routines ;;;;;;;;;;;;;;;;;;;;;;



;;;;;;; file creation routines ;;;;;;;;;;;;;;;;;;;;;;;


undef("create_output_file")
procedure create_output_file(FILES:list,file_path:string,file_name:string,\
									var_array:string)
local FOUT, Times, XLAT, XLONG, XTEMPLATE, dimNames, dimSizes, dimUnlim, \
		var_att, ii, dSizes, dNames, west_east, south_north
begin

;; set NetCDF options, to take advantage of space-saving features
setfileoption("nc","Format","NetCDF4classic")
;setfileoption("nc","CompressionLevel","4")

print("creating file:")
print(file_path+"/"+file_name)
;; open file (this will fail if the file already exists - which is a good thing!)
FOUT  = addfile (file_path+"/"+file_name, "c")
setfileoption(FOUT,"DefineMode",True)

;; copy file attributes from our template files
fileattdef( FOUT, FILES[0] )

;; get the required variables from template files
Times = FILES[0]->Times
XLAT = FILES[0]->XLAT
XLONG = FILES[0]->XLONG
;; save a template variable, to use for creating the variables which need a Time dimensions
XTEMPLATE = XLAT 

;; save only the 1st value for XLAT and XLONG - effectively eliminating the Time dimension from these
;;            variables (this will alter the metadata too for us).
XLAT  := XLAT(0,:,:)
XLONG := XLONG(0,:,:)


;; get the file dimensions, and names
dSizes = getfiledimsizes(FILES[0])
dNames = getvardims(FILES[0])

;; extract west_east and south_north information from the original datafile
west_east   = dSizes(ind(dNames.eq."west_east"))
south_north = dSizes(ind(dNames.eq."south_north"))

;; define the required dimensions for our new file, and add them to the file
dimNames = (/"Time","DateStrLen","west_east","south_north"/)
dimSizes = (/ -1, 19, west_east, south_north /)
dimUnlim = (/True, False, False, False/)
filedimdef(FOUT,dimNames,dimSizes,dimUnlim)

;; define the standard variables, and copy the attributes over
filevardef(FOUT,"Times",typeof(Times),getvardims(Times))
filevardef(FOUT,"XLAT",typeof(XLAT),getvardims(XLAT))
filevardef(FOUT,"XLONG",typeof(XLONG),getvardims(XLONG))
filevarattdef(FOUT,"Times",Times)
filevarattdef(FOUT,"XLAT",XLAT)
filevarattdef(FOUT,"XLONG",XLONG)

;; define the extra variables - for the moment we will assume a lot about these,
;;    such as them all being floats, and 3 dimensional like XLAT (time, w-e, s-n)
var_type = new(dimsizes(var_array),string)
var_type = "float"
filevardef(FOUT,var_array,var_type,getvardims(XTEMPLATE))

;; set up the attributes for the extra variables 
var_att = True
var_att@FieldType   = 104
var_att@MemoryOrder = "XY "
var_att@stagger     = ""
var_att@coordinates = "XLONG XLAT XTIME"
do ii = 0, dimsizes(var_array)-1
	var_att@description = variable_description@$var_array(ii)$
	var_att@units = variable_unit_description@$var_array(ii)$
	filevarattdef(FOUT,var_array(ii),var_att)
end do


;; exit define mode, for writing variable data
setfileoption(FOUT,"DefineMode",False)

;; write the standard data
FOUT->Times = Times
FOUT->XLAT = XLAT
FOUT->XLONG = XLONG


delete(FOUT)

end





end ;; end of script
