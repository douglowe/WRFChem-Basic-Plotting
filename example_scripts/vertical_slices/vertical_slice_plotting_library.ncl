;
; library file for the PLOT_emissions.ncl script
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

load "/Volumes/Odin/YRD_China/data_plotting/scripts/library/variable_loading_library.ncl"

begin

;;;;;;; Setting up the logical variables which we will use to carry 
;;;;;;;   attributes controlling the variable setups
;; used for carrying colour scale information, e.g.:
;;   colour_levels@no2       = (/0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0/) 
colour_levels = True

;; used for carrying information about changing variable units, e.g.:
;;     variable_units@T2  = (/"K","C"/)
variable_units = True


;;; define if the species is a VBS compound or not
VBS_SPECIES = True
VBS_SPECIES@OM_pm1  = True
VBS_SPECIES@OM_pm25 = True
VBS_SPECIES@OM_tot  = True


;;; default colormap - this could be changed in the calling script if needed
colormap = "WhBlGrYeRe"


;;;;;;; settings loading routine
undef("load_var_settings")
procedure load_var_settings(var:string,scenario:string,\
					resources:logical,fixed_colour_scale:logical) 
begin
	
	;;; delete the previous resources
	delete_VarAtts(resources,-1)
	
	;;; set general plotting resources
	resources@MainTitle 	= ""
	resources@InitTime 		= False
	resources@Footer 		= False
	resources@cnFillOn 		= True
	;resources@cnFillMode	= "RasterFill"
	
	if(fixed_colour_scale .and. isatt(colour_levels,var))
		resources@cnLevelSelectionMode = "ExplicitLevels"
		resources@cnLevels = colour_levels@$var$
	end if

end







;;;;;;; variable loading routines
undef("load_emissions")
function load_emissions(FILES:list,var:string)
local variable
begin

	variable = FILES[:]->$var$



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
local var1, var2, var3, var4, var5, var6, vara, dima, itime, var_temp
begin

	;; process the accumulated precipitation data to give hourly precip rates
	if(isStrSubset(var,"_hour"))
		if(var .eq. "RAINC_hour")
			vara = load_variable(FILES,"RAINC",False,variable_units)	
		else if(var .eq. "RAINSH_hour")
			vara = load_variable(FILES,"RAINSH",False,variable_units)	
		else if(var .eq. "RAINNC_hour")
			vara = load_variable(FILES,"RAINNC",False,variable_units)	
		else if(var .eq. "SNOWNC_hour")
			vara = load_variable(FILES,"SNOWNC",False,variable_units)	
		else if(var .eq. "GRAUPELNC_hour")
			vara = load_variable(FILES,"GRAUPELNC",False,variable_units)	
		else if(var .eq. "HAILNC_hour")
			vara = load_variable(FILES,"HAILNC",False,variable_units)	
		else if(var .eq. "rain_hour" .or. var .eq. "precip_hour")
			var1 = load_variable(FILES,"RAINC",False,variable_units)	
			var2 = load_variable(FILES,"RAINSH",False,variable_units)	
			var3 = load_variable(FILES,"RAINNC",False,variable_units)	
			
			if(var .eq. "precip_hour")
				var4 = load_variable(FILES,"SNOWNC",False,variable_units)	
				var5 = load_variable(FILES,"GRAUPELNC",False,variable_units)	
				var6 = load_variable(FILES,"HAILNC",False,variable_units)	
				vara = var1+var2+var3+var4+var5+var6
			else	
				vara = var1+var2+var3
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
	
		;;; call routine for extracting whole domain volume, assume it's not a VBS species
		var_temp = load_variable(FILES,var,False,variable_units)
	
	end if



	return(var_temp)


end





;;;;;;; plotting routines ;;;;;;;;;;;;;;;;;;;;;;;


;;;; map plotting ;;;;;
undef("plot_2d_vertical_slices_from_4d_data")
procedure plot_2d_vertical_slices_from_4d_data(FILES:list,plot_var:numeric,z:numeric,\
		var_name:string,output_directory:string,scenario:string,\
		var_resources:logical,plot_resources:logical,\
		plot_type:string,pivot_point:string,zmax:float)
local times, file_times, xlat, xlon, vardims, loc, angles, x_sect, \
			zmin, nz, itime, idir, var_plane, X_plane, X_desc, \
			dimsX, xmin, xmax, xspan, nx, \
			zz, b, zmax_pos, zspan, \
			output_file_name, output_file, wks, \
			opts_xy, contour, plot
			
begin

	;;; get strings of the time information
	times = wrf_user_getvar(FILES,"times",-1)

	;;; get strings of the time information formatted to make nice file names
	file_times = create_sensible_time_strings(FILES)

	;;; get the latitude and longitude data
	xlat = wrf_user_getvar(FILES,"XLAT",0)
	xlon = wrf_user_getvar(FILES,"XLONG",0)

	;;; get dimensions of variable to plot
	vardims = dimsizes(plot_var)

	;;; pull out the grid points for the cross section
	loc = wrf_user_ll_to_ij(FILES[0], pivot_point@Lon, pivot_point@Lat, True)
	; convert to NCL subscripts
	loc = loc - 1

	;;; define the angle for the cross-section (0 is S-N, 90 is W-E)
	angles = (/ 0.0 , 90.0 /)

	;;; name the cross_sections
	x_sect = (/"S-N_section","W-E_section"/)
	
	;;; define the height labels
	zmin = 0.
	nz   = floattoint(zmax*2 + 1) 

	;;; loop through model times
	do itime = 0,(vardims(0)-1)

		do idir = 0,1
		
			;;; extract the vertical slices
			var_plane = wrf_user_intrp3d(plot_var(itime,:,:,:),z(itime,:,:,:),"v",loc,angles(idir),False)

			;;; get the X-axis info
			if(angles(idir) .eq. 0.0)
				X_plane = wrf_user_intrp2d(xlat,loc,angles(idir),False)
				X_desc = "latitude"
			else if(angles(idir) .eq. 90.0)
				X_plane = wrf_user_intrp2d(xlon,loc,angles(idir),False)			
				X_desc = "longitude"
			end if
			end if
			
			;; X-axis labels
			dimsX = dimsizes(X_plane)
			xmin  = X_plane(0)
			xmax  = X_plane(dimsX(0)-1)
			xspan = dimsX(0)-1
			nx    = floattoint( (xmax-xmin)/2 + 1)


			;;; work out how to trim the data set for our reduced vertical height
	        zz = wrf_user_intrp3d(z(itime,:,:,:),z(itime,:,:,:),"v",loc,angles(idir),False)
    		b = ind(zz(:,0) .gt. zmax*1000. )
          	zmax_pos = b(0) - 1
          	if ( abs(zz(zmax_pos,0)-zmax*1000.) .lt. abs(zz(zmax_pos+1,0)-zmax*1000.) ) then
            	zspan = b(0) - 1
          	else
            	zspan = b(0)
          	end if
          	delete(zz)
          	delete(b)


			;;; name the output file
			output_file_name = var_name+"_"+scenario+"_"+pivot_point+"_"+x_sect(idir)+"_"+file_times(itime)
			output_file = output_directory+"/"+output_file_name

			;;; open plotting file
			wks = gsn_open_wks(plot_type,output_file)       ; Create a plot workstation
	
			;;; define our colourscale
			gsn_define_colormap(wks,colormap)


			opts_xy   = var_resources
			opts_xy@TimeLabel				= times(itime)
			opts_xy@tiXAxisString           = X_desc
			opts_xy@tiYAxisString           = "Height (km)"
			opts_xy@cnMissingValPerimOn     = True
			opts_xy@cnMissingValFillColor   = 0
			opts_xy@cnMissingValFillPattern = 11
			opts_xy@tmYLMode                = "Explicit"
			opts_xy@tmYLValues              = fspan(0,zspan,nz)                    ; Create tick marks
			opts_xy@tmYLLabels              = sprintf("%.1f",fspan(zmin,zmax,nz))  ; Create labels
			opts_xy@tmXBMode                = "Explicit"
			opts_xy@tmXBValues              = fspan(0,xspan,nx)                    ; Create tick marks
			opts_xy@tmXBLabels              = sprintf("%.1f",fspan(xmin,xmax,nx))  ; Create labels
			opts_xy@tmXBLabelFontHeightF    = 0.015
			opts_xy@tiXAxisFontHeightF      = 0.020
			opts_xy@tiYAxisFontHeightF      = 0.020
			opts_xy@tmXBMajorLengthF        = 0.02
			opts_xy@tmYLMajorLengthF        = 0.02
			opts_xy@tmYLLabelFontHeightF    = 0.015
			opts_xy@PlotOrientation         = var_plane@Orientation
 			
			opts_xy@FieldTitle = var_name
			opts_xy@SubFieldTitle = "at UTC time "+times(itime)

			contour = wrf_contour(FILES[itime],wks,var_plane(0:zmax_pos,:),opts_xy)
			plot = wrf_overlays(FILES[itime],wks,contour,plot_resources)
		

			;;; close the plotting file
			delete(wks)		

			;;; delete temporary data array
			delete(var_plane)
			delete(X_plane)
			delete(nx)
			delete(opts_xy)
		
		end do
		
	
	end do

end





end ;; end of script