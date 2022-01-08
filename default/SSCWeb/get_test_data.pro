
FUNCTION Get_Test_Data, satellite1, satellite2, staTime, endTime
; satellite 1
  sat1 = Obj_New( "IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_SATELLITESPECIFICATION", $
    "gov.nasa.gsfc.spdf.ssc.client.SatelliteSpecification" )
  sat1 -> setID, satellite1
  sat1 -> setResolutionFactor, 3
; satellite 2  
  sat2 = Obj_New( "IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_SATELLITESPECIFICATION", $
    "gov.nasa.gsfc.spdf.ssc.client.SatelliteSpecification" )
  sat2 -> setID, satellite2
  sat1 -> setResolutionFactor, 3
; Space region filter options
  spaceRegionsFilterOptions = Obj_New( "IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_SPACEREGIONSFILTEROPTIONS", $
    "gov.nasa.gsfc.spdf.ssc.client.SpaceRegionsFilterOptions" )
  spaceRegionFilterOptions -> setDaysideMagnetosheath, 1
; Hemisphere options  
  hemisphereOptions = Obj_New( "IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_HEMISPHEREOPTIONS", $
    "gov.nasa.gsfc.spdf.ssc.client.HemisphereOptions" )
  hemisphereOptions -> setNorth, 1
  hemisphereOptions -> setSouth, 1
; Mapped region filter options
  mappedRegionFilterOptions = Obj_New( "IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_MAPPEDREGIONFILTEROPTIONS", $
    "gov.nasa.gsfc.spdf.ssc.client.MappedRegionFilterOptions" )
; Region filter options
  regionFilterOptions = Obj_New( "IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_REGIONFILTEROPTIONS", $
    "gov.nasa.gsfc.spdf.ssc.client.RegionFilterOptions" )
  

END