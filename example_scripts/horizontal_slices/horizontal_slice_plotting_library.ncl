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



;;; define if the species is a VBS compound or not
VBS_SPECIES = True
VBS_SPECIES@OM_pm1  = True
VBS_SPECIES@OM_pm25 = True
VBS_SPECIES@OM_tot  = True
VBS_SPECIES@OM_tot_belowcloudbase = True


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



;;;;;;; plotting routines ;;;;;;;;;;;;;;;;;;;;;;;


;;;; map plotting ;;;;;
undef("plot_2d_horizontal_slice_from_4d_data")
procedure plot_2d_horizontal_slice_from_4d_data(FILES:list,plot_var:numeric,\
		var_name:string,output_directory:string,scenario:string,\
		var_resources:logical,map_resources:logical,plot_resources:logical,\
		wind_resources:logical,wind_switch:logical,\
		plot_type:string,altitude_levels:numeric)
local output_file_name, output_file, wks, vardims, times, file_times, itime, ilev, \
		ua, va, alt_km, ua_var, va_var
begin

	;;; get strings of the time information
	times = wrf_user_getvar(FILES,"times",-1)

	;;; get strings of the time information formatted to make nice file names
	file_times = create_sensible_time_strings(FILES)


	;;; get dimensions of variable to plot
	vardims = dimsizes(plot_var)

	if(wind_switch)
		;;; load wind data
		ua = wrf_user_getvar(FILES,"ua",-1)
		va = wrf_user_getvar(FILES,"va",-1)
		
		;;; convert from metres to km for the altitude
		alt_km = altitude_levels / 1000.0


		;;; extract horizontal slices from the domain volume
		ua_var = wrf_user_vert_interp(FILES,ua,"ght_msl",alt_km,True)
		va_var = wrf_user_vert_interp(FILES,va,"ght_msl",alt_km,True)

		;;; replace any altitudes of zero with the data from the ground level
		do ilev = 0,dimsizes(alt_km)-1
			if(alt_km(ilev) .eq. 0.0)
				ua_var(:,ilev,:,:) = (/ua(:,0,:,:)/)
				va_var(:,ilev,:,:) = (/va(:,0,:,:)/)
			end if
		end do

	end if


	;;; loop through model times
	do itime = 0,(vardims(0)-1)

		;;; name the output file
		output_file_name = var_name+"_"+scenario+"_"+file_times(itime)
		output_file = output_directory+"/"+output_file_name

		;;; open plotting file
		wks = gsn_open_wks(plot_type,output_file)       ; Create a plot workstation
	
		;;; define our colourscale
		gsn_define_colormap(wks,colormap)


		;;; plot the model data
		do ilev = 0,dimsizes(altitude_levels)-1
			if(isatt(variable_description,var_name))
				var_resources@FieldTitle = variable_description@$var_name$
			else if(isatt(plot_var,"description"))
				var_resources@FieldTitle = plot_var@description
			else
				var_resources@FieldTitle = var_name
			end if
			end if
			var_resources@SubFieldTitle = "at UTC time "+times(itime)
			if(altitude_levels(ilev) .eq. 0.0)
				var_resources@PlotLevelID = "ground level"
			else
				var_resources@PlotLevelID = flt2string(altitude_levels(ilev))+" metres ASL"
			end if
			
			; check to make sure there's *some* data - if not set all to zero
			if( all( ismissing(plot_var(itime,ilev,:,:)) ) )
				print("no data for level "+ilev+" --- setting all to zero for this level")
				plot_var(itime,ilev,:,:) = 0.0
			end if
			
			contour = wrf_contour(FILES[0],wks,plot_var(itime,ilev,:,:),var_resources)

			if(wind_switch) ; plot wind vectors, and combine with contour plot
				vector = wrf_vector(FILES[0],wks,ua_var(itime,ilev,:,:),va_var(itime,ilev,:,:),wind_resources)		
				plot = wrf_map_overlays(FILES[0],wks,(/contour,vector/),plot_resources,map_resources)
			else ; only plot the contours
				plot = wrf_map_overlays(FILES[0],wks,(/contour/),plot_resources,map_resources)
			end if
		end do
	
		;;; close the plotting file
		delete(wks)		
	
	end do

end


;;;; map plotting ;;;;;
undef("plot_2d_data")
procedure plot_2d_data(FILES:list,plot_var:numeric,\
		var_name:string,output_directory:string,scenario:string,\
		var_resources:logical,map_resources:logical,plot_resources:logical,\
		wind_resources:logical,wind_switch:logical,plot_type:string)
local output_file_name, output_file, wks, vardims, times, file_times, itime, ilev
begin

	;;; get strings of the time information
	times = wrf_user_getvar(FILES,"times",-1)
	
	;;; get strings of the time information formatted to make nice file names
	file_times = create_sensible_time_strings(FILES)

	;;; get dimensions of variable to plot
	vardims = dimsizes(plot_var)

	if(wind_switch)
		;;; load wind data
		u10 = wrf_user_getvar(FILES,"U10",-1)
		v10 = wrf_user_getvar(FILES,"V10",-1)
	end if


	;;; loop through model times
	do itime = 0,(vardims(0)-1)

		;;; name the output file
		output_file_name = var_name+"_"+scenario+"_"+file_times(itime)
		output_file = output_directory+"/"+output_file_name

		;;; open plotting file
		wks = gsn_open_wks(plot_type,output_file)       ; Create a plot workstation
	
		;;; define our colourscale
		gsn_define_colormap(wks,colormap)


		;;; plot the model data
		if(isatt(variable_description,var_name))
			var_resources@FieldTitle = variable_description@$var_name$
		else if(isatt(plot_var,"description"))
			var_resources@FieldTitle = plot_var@description
		else
			var_resources@FieldTitle = var_name
		end if
		end if
		var_resources@SubFieldTitle = "at UTC time "+times(itime)
		;var_resources@PlotLevelID = ""
		
		if( all( ismissing(plot_var(itime,:,:)) ) )
			print("no data --- setting all to zero")
			plot_var(itime,:,:) = 0.0
		end if

		contour = wrf_contour(FILES[0],wks,plot_var(itime,:,:),var_resources)
		
		if(wind_switch) ; plot wind vectors, and combine with contour plot
			vector = wrf_vector(FILES[0],wks,u10(itime,:,:),v10(itime,:,:),wind_resources)		
			plot = wrf_map_overlays(FILES[0],wks,(/contour,vector/),plot_resources,map_resources)
		else ; only plot the contours
			plot = wrf_map_overlays(FILES[0],wks,contour,plot_resources,map_resources)		
		end if
	
		;;; close the plotting file
		delete(wks)		
	
	end do

end




end ;; end of script
