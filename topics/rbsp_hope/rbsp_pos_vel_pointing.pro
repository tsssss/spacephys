;rbsp_efw_position_velocity_crib

;Loads and plots RBSP (Van Allen probes) position and velocity data
;	MLT
;	Mlat
;	Lshell
;	Position (GSE)
;	Velocity (GSE)
;
;Need to have the SPICE ICY software package installed
;
;
;Written by Aaron Breneman, UMN, Dec 2012
;			email: awbrenem@gmail.com





; initialize RBSP environment
	rbsp_efw_init
	!rbsp_efw.user_agent = ''


; set desired probe
	probe = 'a'


; set time of interest to a single day
	date = '2012-10-13'	; UTC.
	duration = 1 ; days.
	timespan, date, duration



	type='calibrated' ; use ADC (raw) numbers or physical (calibrated) units?


; Load definitive spice data
	rbsp_load_spice_kernels



	;Load state data
	rbsp_load_spice_state,probe='a',coord='gse',/no_spice_load  
	rbsp_load_spice_state,probe='b',coord='gse',/no_spice_load  
	
	get_data,'rbspa_state_pos_gse',data=pos_gse_a
	get_data,'rbspb_state_pos_gse',data=pos_gse_b


	;Create position and velocity magnitude variables
	get_data,'rbspa_state_vel_gse',data=vela
	get_data,'rbspb_state_vel_gse',data=velb
	
	vmaga = sqrt(vela.y[*,0]^2 + vela.y[*,1]^2 + vela.y[*,2]^2)
	vmagb = sqrt(velb.y[*,0]^2 + velb.y[*,1]^2 + velb.y[*,2]^2)
	
	store_data,'vmag_a',data={x:vela.x,y:vmaga}
	store_data,'vmag_b',data={x:velb.x,y:vmagb}
	

	if probe eq 'a' then pos_gse = pos_gse_a else pos_gse = pos_gse_b
	rad = sqrt(pos_gse.y[*,0]^2 + pos_gse.y[*,1]^2 + pos_gse.y[*,2]^2)/6370.

	rad_a = sqrt(pos_gse_a.y[*,0]^2 + pos_gse_a.y[*,1]^2 + pos_gse_a.y[*,2]^2)/6370.
	rad_b = sqrt(pos_gse_b.y[*,0]^2 + pos_gse_b.y[*,1]^2 + pos_gse_b.y[*,2]^2)/6370.


	store_data,'radius',data={x:pos_gse.x,y:rad}
	options,'rbsp'+probe+'_state_pos_gse','panel_size',1
	options,'rbsp'+probe+'_state_vel_gse','panel_size',1
	options,'radius','panel_size',0.6
	options,'radius','ytitle','Rad!C[RE]'
	


	cotrans,'rbspa_state_pos_gse','rbspa_state_pos_gsm',/GSE2GSM	
	cotrans,'rbspb_state_pos_gse','rbspb_state_pos_gsm',/GSE2GSM	

	cotrans,'rbspa_state_vel_gse','rbspa_state_vel_gsm',/GSE2GSM	
	cotrans,'rbspb_state_vel_gse','rbspb_state_vel_gsm',/GSE2GSM	


;For calculating Mlat only
	cotrans,'rbspa_state_pos_gsm','rbspa_state_pos_sm',/GSM2SM	
	cotrans,'rbspb_state_pos_gsm','rbspb_state_pos_sm',/GSM2SM	







;Calculate magnetic latitude	

	get_data,'rbspa_state_pos_sm',data=pos_sm_a
	get_data,'rbspb_state_pos_sm',data=pos_sm_b

	
	dr2a = sqrt(pos_sm_a.y[*,0]^2 + pos_sm_a.y[*,1]^2)
	dz2a = pos_sm_a.y[*,2]
	dr2b = sqrt(pos_sm_b.y[*,0]^2 + pos_sm_b.y[*,1]^2)
	dz2b = pos_sm_b.y[*,2]

	mlat_a = atan(dz2a,dr2a)
	mlat_b = atan(dz2b,dr2b)
	
	if probe eq 'a' then mlat = mlat_a else mlat = mlat_b


	store_data,'mlat',data={x:pos_gse.x,y:mlat/!dtor}
	options,'mlat','format','(f5.1)'        
	options,'mlat','ytitle','mlat'






;Calculate L-shell



	;Method 1
	Lshell_a = rad_a/(cos(!dtor*mlat_a)^2)  ;L-shell in centered dipole
	Lshell_b = rad_b/(cos(!dtor*mlat_b)^2)  ;L-shell in centered dipole


	if probe eq 'a' then Lshell = Lshell_a else Lshell = Lshell_b



	;Method 2
	;Position data must be in km
	ttrace2equator,'rbspa_state_pos_gsm',newname='rbspa_out_foot',/km
	ttrace2equator,'rbspb_state_pos_gsm',newname='rbspb_out_foot',/km
	

	get_data,'rbspa_out_foot',data=d
	Lshell_a = sqrt(d.y[*,0]^2 + d.y[*,1]^2 + d.y[*,2]^2)/6370.
	get_data,'rbspb_out_foot',data=d
	Lshell_b = sqrt(d.y[*,0]^2 + d.y[*,1]^2 + d.y[*,2]^2)/6370.

	if probe eq 'a' then Lshell = Lshell_a else Lshell = Lshell_b


	store_data,'lshell',data={x:pos_gse.x,y:lshell}
	options,'lshell','format','(f5.1)'
	options,'lshell','ytitle','Lshell'



	







	

;Calculate MLT
	
	angle_tmp = atan(pos_gse_a.y[*,1],pos_gse_a.y[*,0])/!dtor
	goo = where(angle_tmp lt 0.)
	if goo[0] ne -1 then angle_tmp[goo] = 360. - abs(angle_tmp[goo])
	angle_rad_a = angle_tmp * 12/180. + 12.
	goo = where(angle_rad_a ge 24.)
	if goo[0] ne -1 then angle_rad_a[goo] = angle_rad_a[goo] - 24

	angle_tmp = atan(pos_gse_b.y[*,1],pos_gse_b.y[*,0])/!dtor
	goo = where(angle_tmp lt 0.)
	if goo[0] ne -1 then angle_tmp[goo] = 360. - abs(angle_tmp[goo])
	angle_rad_b = angle_tmp * 12/180. + 12.
	goo = where(angle_rad_b ge 24.)
	if goo[0] ne -1 then angle_rad_b[goo] = angle_rad_b[goo] - 24

	if probe eq 'a' then angle_rad = angle_rad_a else angle_rad = angle_rad_b
	
	
	store_data,'mlt',data={x:pos_gse.x,y:angle_rad}




	



;Find differences in MLT and L b/t the two sc
	store_data,'mlt_a',data={x:pos_gse_a.x,y:angle_rad_a}
	store_data,'lshell_a',data={x:pos_gse_a.x,y:lshell_a}
	store_data,'mlat_a',data={x:pos_gse_a.x,y:mlat_a/!dtor}

	store_data,'mlt_b',data={x:pos_gse_b.x,y:angle_rad_b}
	store_data,'lshell_b',data={x:pos_gse_b.x,y:lshell_b}
	store_data,'mlat_b',data={x:pos_gse_b.x,y:mlat_b/!dtor}
	


	;Interpolate to get GSE pos of both sc on the same times.
	tinterpol_mxn,'mlt_a',pos_gse.x,newname='mlt_a2'
	tinterpol_mxn,'mlt_b',pos_gse.x,newname='mlt_b2'	
	dif_data,'mlt_a2','mlt_b2',newname='mlt_diff'

	tinterpol_mxn,'lshell_a',pos_gse.x,newname='lshell_a2'
	tinterpol_mxn,'lshell_b',pos_gse.x,newname='lshell_b2'	
	dif_data,'lshell_a2','lshell_b2',newname='lshell_diff'

	tinterpol_mxn,'mlat_a',pos_gse.x,newname='mlat_a2'
	tinterpol_mxn,'mlat_b',pos_gse.x,newname='mlat_b2'	
	dif_data,'mlat_a2','mlat_b2',newname='mlat_diff'


	options,'mlat_a2','ytitle','MLATa!Chours'
	options,'mlat_b2','ytitle','MLATb!Chours'
	options,'mlt_diff','ytitle','MLTa-!CMLTb!Chours'
	options,'lshell_diff','ytitle','LSHELLa-!CLSHELLb'
	options,'mlat_diff','ytitle','MLATa-!CMLATb!Cdegrees'
	options,'rbspa_state_vel_gse','ytitle','RBSPa!Cvelocity!CGSE'
	options,'rbspb_state_vel_gse','ytitle','RBSPb!Cvelocity!CGSE'


	ylim,'mlat_diff',-10,10
	ylim,'mlt_diff',-5,5

	options,'vmag_a','ytitle','|V| RBSPa!Ckm/s'
	options,'vmag_b','ytitle','|V| RBSPb!Ckm/s'





;Create x-axis labels for the plots
	var_label = ['lshell','mlat','mlt']





;Plot various quantities

	tplot_options,'title','rbsp position data - ' + date		


	;Plot position quantities
	tplot,['rbspa_state_pos_gse',$
		   'rbspb_state_pos_gse',$
		   'rbspa_state_pos_gsm',$
		   'rbspb_state_pos_gsm'],$ 
			var_label = var_label


	tplot,['mlat_a2',$
		   'mlat_b2',$
		   'lshell_a2',$
		   'lshell_b2',$
		   'mlt_a2',$
		   'mlt_b2'],$ 
			var_label = var_label



	
	;Plot velocity quantities
	tplot,['vmag_a',$
		   'rbspa_state_vel_gse',$
		   'rbspa_state_vel_gsm',$
		   'vmag_b',$
		   'rbspb_state_vel_gse',$
		   'rbspb_state_vel_gsm'],$ 
			var_label = var_label





;--------------------------------------------------------------------
;Get the Wgse direction for the MGSE transformation
;--------------------------------------------------------------------

	
		time2=time_double(date) ; first get unix time double for beginning of day
	
		;Grab first and last time of day
		time2 = [time2,time2+86399.]
		time3=time_string(time2, prec=6) ; turn it back into a string for ISO conversion
		strput,time3,'T',10 ; convert TPLOT time string 'yyyy-mm-dd/hh:mm:ss.msec' to ISO 'yyyy-mm-ddThh:mm:ss.msec'
		cspice_str2et,time3,et2 ; convert ISO time string to SPICE ET
	
		
		cspice_pxform,'RBSP'+strupcase(probe)+'_SCIENCE','GSE',et2[0],pxform1
		cspice_pxform,'RBSP'+strupcase(probe)+'_SCIENCE','GSE',et2[1],pxform2
		
		wsc1=dblarr(3)
		wsc1[2]=1d
		wsc_GSE1=dblarr(3)
		
		wsc2=dblarr(3)
		wsc2[2]=1d
		wsc_GSE2=dblarr(3)
		
		; Calculate the modified MGSE directions	
		wsc_GSE1 = pxform1 ## wsc1  ;start of day
		wsc_GSE2 = pxform2 ## wsc2  ;end of day	
	
	


;--------------------------------------------------------------------
;Use the Wgse direction to rotate from GSE to MGSE
;--------------------------------------------------------------------

	
	
	
		;Puts velocity data in MGSE
;		rbsp_gse2mgse,rbspx+'_state_vel_gse',reform(wsc_GSE1),newname=rbspx+'_vel_mgse'
		rbsp_gse2mgse,rbspx+'_vel_gse',reform(wsc_GSE1),newname=rbspx+'_vel_mgse'


		;Puts EMFISIS L3 data in MGSE
		get_data,rbspx+'_mag_gse_from_l3',data=mag_gse
		if is_struct(mag_gse) then $
			rbsp_gse2mgse,rbspx+'_mag_gse_from_l3',reform(wsc_GSE1),newname=rbspx+'_mag_mgse_from_l3'


		options,rbspx+'_mag_mgse_from_l3','ytitle','Bfield!CMGSE'
		options,rbspx+'_mag_mgse_from_l3','ysubtitle','[nT]'
		options,rbspx+'_mag_mgse_from_l3','colors',[2,4,6]

		options,rbspx+'_mag_mgse_from_uvw','ytitle','Bfield!CMGSE'
		options,rbspx+'_mag_mgse_from_uvw','ysubtitle','[nT]'
		options,rbspx+'_mag_mgse_from_uvw','colors',[2,4,6]
			
        endif

end