;
; NOSA HEADER START
;
; The contents of this file are subject to the terms of the NASA Open 
; Source Agreement (NOSA), Version 1.3 only (the "Agreement").  You may 
; not use this file except in compliance with the Agreement.
;
; You can obtain a copy of the agreement at
;   docs/NASA_Open_Source_Agreement_1.3.txt
; or 
;   http://sscweb.gsfc.nasa.gov/WebServices/NASA_Open_Source_Agreement_1.3.txt.
;
; See the Agreement for the specific language governing permissions
; and limitations under the Agreement.
;
; When distributing Covered Code, include this NOSA HEADER in each
; file and include the Agreement file at 
; docs/NASA_Open_Source_Agreement_1.3.txt.  If applicable, add the 
; following below this NOSA HEADER, with the fields enclosed by 
; brackets "[]" replaced with your own identifying information: 
; Portions Copyright [yyyy] [name of copyright owner]
;
; NOSA HEADER END
;
; Copyright (c) 2007-2008 United States Government as represented by the 
; National Aeronautics and Space Administration. No copyright is claimed 
; in the United States under Title 17, U.S.Code. All Other Rights Reserved.
;
; $Id: WsExample.pro,v 1.2.4.2.2.4 2008/06/17 17:27:50 bharris Exp $
;


FUNCTION getTestOutputOptions

    ; get the value of the Java enumeration for the coordinate system
    ; we want in an IDL variable
    ;
    CoordinateSystem = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_COORDINATESYSTEM', $
        'gov.nasa.gsfc.spdf.ssc.client.CoordinateSystem')
    CoordinateSystem -> getProperty, GEO=geoEnum

    ; get the values of the Java enumeration for the coordinate components
    ; we want in IDL variables
    ;
    CoordinateComponent = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_COORDINATECOMPONENT', $
        'gov.nasa.gsfc.spdf.ssc.client.CoordinateComponent')
    CoordinateComponent -> getProperty, X=xComponentEnum
    CoordinateComponent -> getProperty, Y=yComponentEnum
    CoordinateComponent -> getProperty, Z=zComponentEnum
    CoordinateComponent -> getProperty, LAT=latComponentEnum
    CoordinateComponent -> getProperty, LON=lonComponentEnum
    CoordinateComponent -> getProperty, LOCAL_TIME=localTimeComponentEnum

    ; create the FilteredCoordinateOptions objects we need to represent
    ; coordinate options we want
    ;
    geoXOption = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_FILTEREDCOORDINATEOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.FilteredCoordinateOptions')
    geoXOption -> setCoordinateSystem, geoEnum
    geoXOption -> setComponent, xComponentEnum
    geoXOption -> setFilter, OBJ_NEW()

    geoYOption = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_FILTEREDCOORDINATEOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.FilteredCoordinateOptions')
    geoYOption -> setCoordinateSystem, geoEnum
    geoYOption -> setComponent, yComponentEnum
    geoYOption -> setFilter, OBJ_NEW()

    geoZOption = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_FILTEREDCOORDINATEOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.FilteredCoordinateOptions')
    geoZOption -> setCoordinateSystem, geoEnum
    geoZOption -> setComponent, zComponentEnum
    geoZOption -> setFilter, OBJ_NEW()

    geoLatOption = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_FILTEREDCOORDINATEOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.FilteredCoordinateOptions')
    geoLatOption -> setCoordinateSystem, geoEnum
    geoLatOption -> setComponent, latComponentEnum
    geoLatOption -> setFilter, OBJ_NEW()

    geoLonOption = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_FILTEREDCOORDINATEOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.FilteredCoordinateOptions')
    geoLonOption -> setCoordinateSystem, geoEnum
    geoLonOption -> setComponent, lonComponentEnum
    geoLonOption -> setFilter, OBJ_NEW()

    geoLocalTimeOption = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_FILTEREDCOORDINATEOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.FilteredCoordinateOptions')
    geoLocalTimeOption -> setCoordinateSystem, geoEnum
    geoLocalTimeOption -> setComponent, localTimeComponentEnum
    geoLocalTimeOption -> setFilter, OBJ_NEW()


    ; get the Java Hemisphere enumeration values we need in IDL variables
    ;
    HemisphereEnum = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_HEMISPHERE', $
        'gov.nasa.gsfc.spdf.ssc.client.Hemisphere')
    HemisphereEnum -> getProperty, NORTH=northHemisphereEnum
    HemisphereEnum -> getProperty, SOUTH=southHemisphereEnum

    ; create the BFieldTraceOptions we want
    ;
    geoNorthBFieldTrace = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_BFIELDTRACEOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.BFieldTraceOptions')

    geoNorthBFieldTrace -> setCoordinateSystem, geoEnum
    geoNorthBFieldTrace -> setFieldLineLength, 1
    geoNorthBFieldTrace -> setFootpointLatitude, 1
    geoNorthBFieldTrace -> setFootpointLongitude, 1
    geoNorthBFieldTrace -> setHemisphere, northHemisphereEnum

;    geoSouthBFieldTrace = OBJ_NEW( $
;        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_BFIELDTRACEOPTIONS', $
;        'gov.nasa.gsfc.spdf.ssc.client.BFieldTraceOptions')
;
;    geoSouthBFieldTrace -> setCoordinateSystem, geoEnum
;    geoSouthBFieldTrace -> setFieldLineLength, 1
;    geoSouthBFieldTrace -> setFootpointLatitude, 1
;    geoSouthBFieldTrace -> setFootpointLongitude, 1
;    geoSouthBFieldTrace -> setHemisphere, southHemisphereEnum

    ; create the OutputOptions we want
    ;
    outputOptions = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_OUTPUTOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.OutputOptions')

    ; set the OutputOptions
    ;
    outputOptions -> setAllLocationFilters, 1

    coordOptions = outputOptions -> getCoordinateOptions()
    coordOptions -> add, geoXOption
    coordOptions -> add, geoYOption
    coordOptions -> add, geoZOption
    coordOptions -> add, geoLatOption
    coordOptions -> add, geoLonOption
    coordOptions -> add, geoLocalTimeOption

    outputOptions -> setMinMaxPoints, 2

    bFieldTraceOptions = outputOptions -> getBFieldTraceOptions()

    bFieldTraceOptions -> add, geoNorthBFieldTrace
;    bFieldTraceOptions -> add, geoSouthBFieldTrace

RETURN, outputOptions

END


FUNCTION getTestFormatOptions

    ; create IDL variables representing the Java constants we need
    ;
    dateFormat = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_DATEFORMAT', $
        'gov.nasa.gsfc.spdf.ssc.client.DateFormat')
    dateFormat -> getProperty, YYYY_DDD=yyyydddEnum

    timeFormat = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_TIMEFORMAT', $
        'gov.nasa.gsfc.spdf.ssc.client.TimeFormat')
    timeFormat -> getProperty, HH_MM=hhmmEnum

    distanceUnits = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_DISTANCEUNITS', $
        'gov.nasa.gsfc.spdf.ssc.client.DistanceUnits')
    distanceUnits -> getProperty, RE=reEnum

    degreeFormat = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_DEGREEFORMAT', $
        'gov.nasa.gsfc.spdf.ssc.client.DegreeFormat')
    degreeFormat -> getProperty, DECIMAL=decimalEnum

    latLonFormat = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_LATLONFORMAT', $
        'gov.nasa.gsfc.spdf.ssc.client.LatLonFormat')
    latLonFormat -> getProperty, LAT_90_LON_360=lat90lon360Enum

    ;
    ; now create a FormatOptions object and set its properties using
    ; the values obtained above
    ;
    formatOptions = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_FORMATOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.FormatOptions')

    formatOptions -> setDateFormat, yyyydddEnum
    formatOptions -> setTimeFormat, hhmmEnum
    formatOptions -> setDistanceUnits, reEnum
    formatOptions -> setDistancePrecision, 2
    formatOptions -> setDegreeFormat, decimalEnum
    formatOptions -> setDegreePrecision, 2
    formatOptions -> setLatLonFormat, lat90lon360Enum
    formatOptions -> setCdf, 0
    formatOptions -> setLinesPerPage, 1

RETURN, formatOptions

END


FUNCTION getTestDataRequest, startTime, endTime

    fastSat = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_SATELLITESPECIFICATION', $
        'gov.nasa.gsfc.spdf.ssc.client.SatelliteSpecification')
                                       ;// fast satellite spec.
    fastSat -> setId, 'fast'
    fastSat -> setResolutionFactor, 2

    moonSat = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_SATELLITESPECIFICATION', $
        'gov.nasa.gsfc.spdf.ssc.client.SatelliteSpecification')
                                       ;// fast satellite spec.

; The moon is too far away from fast for iplot to make a very good 
; display so use cluster1 instead.
;    moonSat -> setId, 'moon'
    moonSat -> setId, 'cluster1old'
    moonSat -> setResolutionFactor, 2

    spaceRegionsFilter = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_SPACEREGIONSFILTEROPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.SpaceRegionsFilterOptions')

        
    spaceRegionsFilter -> setDaysideMagnetosheath, 1
    spaceRegionsFilter -> setDaysideMagnetosphere, 1
    spaceRegionsFilter -> setDaysidePlasmasphere, 1
    spaceRegionsFilter -> setHighLatitudeBoundaryLayer, 1
    spaceRegionsFilter -> setInterplanetaryMedium, 1
    spaceRegionsFilter -> setLowLatitudeBoundaryLayer, 1
    spaceRegionsFilter -> setNightsideMagnetosheath, 1
    spaceRegionsFilter -> setNightsideMagnetosphere, 1
    spaceRegionsFilter -> setNightsidePlasmasphere, 1
    spaceRegionsFilter -> setPlasmaSheet, 1
    spaceRegionsFilter -> setTailLobe, 1

    hemisphereOptions = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_HEMISPHEREOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.HemisphereOptions')

    hemisphereOptions -> setNorth, 1
    hemisphereOptions -> setSouth, 1

    radialTraceRegionsFilter = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_MAPPEDREGIONFILTEROPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.MappedRegionFilterOptions')
        
    radialTraceRegionsFilter -> setCusp, hemisphereOptions
    radialTraceRegionsFilter -> setCleft, hemisphereOptions
    radialTraceRegionsFilter -> setAuroralOval, hemisphereOptions
    radialTraceRegionsFilter -> setPolarCap, hemisphereOptions
    radialTraceRegionsFilter -> setMidLatitude, hemisphereOptions
    radialTraceRegionsFilter -> setLowLatitude, 1

    magneticTraceRegionsFilter = radialTraceRegionsFilter

    regionFilters = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_REGIONFILTEROPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.RegionFilterOptions')

    regionFilters -> setSpaceRegions, spaceRegionsFilter
    regionFilters -> setRadialTraceRegions, radialTraceRegionsFilter
    regionFilters -> setMagneticTraceRegions, magneticTraceRegionsFilter

    locationFilter = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_LOCATIONFILTER', $
        'gov.nasa.gsfc.spdf.ssc.client.LocationFilter')

    lowerLimit = OBJ_NEW( $
        'IDLJavaObject$JAVA_LANG_DOUBLE', 'java.lang.Double', -500.0D)

    upperLimit = OBJ_NEW( $
        'IDLJavaObject$JAVA_LANG_DOUBLE', 'java.lang.Double', 500.0D)

    locationFilter -> setMinimum, 1
    locationFilter -> setMaximum, 1
    locationFilter -> setLowerLimit, lowerLimit
    locationFilter -> setUpperLimit, upperLimit

    locationFilterOptions = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_LOCATIONFILTEROPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.LocationFilterOptions')

    locationFilterOptions -> setAllFilters, 1
    locationFilterOptions -> setDistanceFromCenterOfEarth, locationFilter
    locationFilterOptions -> setMagneticFieldStrength, locationFilter
    locationFilterOptions -> setDistanceFromNeutralSheet, locationFilter
    locationFilterOptions -> setDistanceFromBowShock, locationFilter
    locationFilterOptions -> setDistanceFromMagnetopause, locationFilter
    locationFilterOptions -> setDipoleLValue, locationFilter
    locationFilterOptions -> setDipoleInvariantLatitude, locationFilter

    externalBFieldModel = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_BFIELDMODELPARAMETERS', $
        'gov.nasa.gsfc.spdf.ssc.client.BFieldModelParameters')

    externalBFieldModelEnum = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_EXTERNALBFIELDMODEL', $
        'gov.nasa.gsfc.spdf.ssc.client.ExternalBFieldModel')
    externalBFieldModelEnum -> getProperty, T_96=t96ExternalBFieldModel

    externalBFieldModel -> setModel, t96ExternalBFieldModel
    externalBFieldModel -> setUseFixedValues, 1
    externalBFieldModel -> setSolarWindPressure, 2.1
    externalBFieldModel -> setDst, -20
    externalBFieldModel -> setByImf, 0.0
    externalBFieldModel -> setBzImf, 0.0

    bFieldModelOptions = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_BFIELDMODELOPTIONS', $
        'gov.nasa.gsfc.spdf.ssc.client.BFieldModelOptions')

    internalBFieldModelEnum = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_INTERNALBFIELDMODEL', $
        'gov.nasa.gsfc.spdf.ssc.client.InternalBFieldModel')
    internalBFieldModelEnum -> getProperty, IGRF=igrfInternalBFieldModel

    bFieldModelOptions -> setInternalModel, igrfInternalBFieldModel
    bFieldModelOptions -> setExternalModel, externalBFieldModel
    bFieldModelOptions -> setFieldLinesStopAltitude, 200.0


    dataRequest = OBJ_NEW( $
        'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_DATAREQUEST', $
        'gov.nasa.gsfc.spdf.ssc.client.DataRequest')

    sats = dataRequest -> getSatellites()
    sats -> add, fastSat
    sats -> add, moonSat

    outputOptions = getTestOutputOptions()

    dataRequest -> setOutputOptions, outputOptions

    dataRequest -> setBeginTime, startTime
    dataRequest -> setEndTime, endTime

;    dataRequest -> setRegionFilterOptions, regionFilters
;    dataRequest -> setLocationFilterOptions, locationFilterOptions
;    dataRequest -> setBFieldModelOptions, bFieldModelOptions

; FormatOptions are only applicable to a DataFileRequest (that produces
; a CDF or text file).
;    testFormatOptions = getTestFormatOptions()
;     dataRequest -> setFormatOptions, testFormatOptions

RETURN, dataRequest

END


FUNCTION kmToRe, value

Re = 6378.160
RETURN, value/Re

END


PRO WsExample

; In many instances, the following error handling causes more trouble 
; (e.g., sending IDL into an infinite loop) than the value of the debugging 
; information it is intended to provide.  But it may be helpful in some
; cases.
;
;oJBridgeSession = OBJ_NEW('IDLJavaObject$IDLJAVABRIDGESESSION')
;CATCH, error_status
;IF (error_status NE 0) THEN BEGIN

;    oJExc = oJBridgeSession -> GetException()
;    HELP, oJExc
;    PRINT, 'Exception thrown:', oJExc -> ToString()
;    oJExc -> PrintStackTrace
;    OBJ_DESTROY, oJExc
;    RETURN
;    CATCH, /CANCEL
;ENDIF

;// ---------- SSC Web Services code ----------
system = OBJ_NEW('IDLJavaObject$Static$JAVA_LANG_SYSTEM', $
                 'java.lang.System')
oldHttpAgent = system -> setProperty('http.agent', 'IdlWsExample2')

wsdlUrl = OBJ_NEW('IDLJavaObject$JAVA_NET_URL', 'java.net.URL', $
    'http://sscweb.gsfc.nasa.gov/WS/ssc/2/SatelliteSituationCenterService?wsdl')

sscWsQName = OBJ_NEW('IDLJavaObject$JAVAX_XML_NAMESPACE_QNAME', $
    'javax.xml.namespace.QName', 'http://ssc.spdf.gsfc.nasa.gov/', $
    'SatelliteSituationCenterService')

satelliteSituationCenterService = OBJ_NEW( $
    'IDLJavaObject$GOV_NASA_GSFC_SPDF_SSC_CLIENT_SATELLITESITUATIONCENTERSERVICE', $
    'gov.nasa.gsfc.spdf.ssc.client.SatelliteSituationCenterService', $
    wsdlUrl, sscWsQName)

ssc = satelliteSituationCenterService -> getSatelliteSituationCenterPort()


sats = ssc -> getAllSatellites()
;// ---------- end SSC Web Services code ----------

PRINT, 'Satellites:'

FOR i=0, sats -> size() - 1 DO BEGIN

  sat = sats -> get(i)

  startTime = sat -> getStartTime()
  endTime = sat -> getEndTime()

  PRINT, '  ', sat -> getId(), '  ', sat -> getName(), $
         '  ', sat -> getResolution(), 's  ', startTime -> toString(), $
         ' - ', endTime -> toString()
; PRINT, '  ', sat -> getGeometry(), '  ', sat -> getTrajectoryGeometry()
ENDFOR

;// ---------- SSC Web Services code ----------
stations = ssc -> getAllGroundStations()
;// ---------- end SSC Web Services code ----------

PRINT, 'Ground Stations:'

FOR i=0, stations -> size() - 1 DO BEGIN

  station = stations -> get(i)

  PRINT, '  ', station -> getId(), '  ', station -> getLatitude(), $
         '  ', station -> getLongitude(), '  ', station -> getName()

ENDFOR

staticDatatypeFactory = OBJ_NEW( $
    'IDLJavaObject$Static$JAVAX_XML_DATATYPE_DATATYPEFACTORY', $
    'javax.xml.datatype.DatatypeFactory')

datatypeFactory = staticDatatypeFactory -> newInstance()

testStart = datatypeFactory -> newXMLGregorianCalendar( $
    2001, 03, 17, 00, 0, 0, 0, 0)

testEnd = datatypeFactory -> newXMLGregorianCalendar( $
    2001, 03, 18, 00, 0, 0, 0, 0)

PRINT, 'testStart = ', testStart -> toString(), ', testEnd = ', $
       testEnd -> toString()

dataRequest = getTestDataRequest(testStart, testEnd)

;// ---------- SSC Web Services code ----------
dataResult = ssc -> getData(dataRequest)
;// ---------- end SSC Web Services code ----------


resultStatusCodeEnum = OBJ_NEW( $
        'IDLJavaObject$Static$GOV_NASA_GSFC_SPDF_SSC_CLIENT_RESULTSTATUSCODE', $
        'gov.nasa.gsfc.spdf.ssc.client.ResultStatusCode')
resultStatusCodeEnum -> getProperty, ERROR=errorStatusCode

resultStatusCode = dataResult -> getStatusCode()

error = resultStatusCode -> equals(errorStatusCode)

IF (error) THEN BEGIN

    MESSAGE, 'ERROR: getData failed'
ENDIF

data = dataResult -> getData()

FOR i=0, data -> size() - 1 DO BEGIN

    ; get i'th satellite's data
    ;
    satData = data -> get(i)

    print, 'Data for ', satData -> getId()

    satTime = satData -> getTime()

    bFieldTraceData = satData -> getBTraceData()

    FOR j=0, bFieldTraceData -> size() - 1 DO BEGIN

        bFieldTrace = bFieldTraceData -> get(j)

        bFieldCoord = bFieldTrace -> getCoordinateSystem()
        bFieldHemisphere = bFieldTrace -> getHemisphere()

        print, 'B Field Trace CoordinateSystem ', bFieldCoord -> toString()
        print, 'B Field Trace Hemisphere ', bFieldHemisphere -> toString()

        latList = bFieldTrace -> getLatitude()
        lonList = bFieldTrace -> getLongitude()

        ; convert the Java Lists of Float objects into primative IDL 
        ; arrays of float values
        ;
        lat = FLTARR(latList -> size(), /NOZERO)
        lon = FLTARR(lonList -> size(), /NOZERO)

        FOR k=0, latList -> size() - 1 DO BEGIN

            latJavaFloat = latList -> get(k)
            lonJavaFloat = lonList -> get(k)

            lat[k] = latJavaFloat -> floatValue()
            lon[k] = lonJavaFloat -> floatValue()

            print, lat[k], lon[k]
        ENDFOR
    ENDFOR

    satCoordData = satData -> getCoordinates()

    satCoordData0 = satCoordData -> get(0)

    satXList = satCoordData0 -> getX()
    satYList = satCoordData0 -> getY()
    satZList = satCoordData0 -> getZ()
    satLatList = satCoordData0 -> getLatitude()
    satLonList = satCoordData0 -> getLongitude()
    satLTList = satCoordData0 -> getLocalTime()

    ; convert the Java Lists of Double and Float objects into primative IDL 
    ; arrays of double and float values
    ;
    x = DBLARR(satXList -> size(), /NOZERO)
    y = DBLARR(satXList -> size(), /NOZERO)
    z = DBLARR(satXList -> size(), /NOZERO)
    lat = FLTARR(satLatList -> size(), /NOZERO)
    lon = FLTARR(satLonList -> size(), /NOZERO)
    localTime = DBLARR(satLTList -> size(), /NOZERO)

    FOR j=0, satXList -> size() - 1 DO BEGIN

        xJavaDouble = satXList -> get(j)
        yJavaDouble = satYList -> get(j)
        zJavaDouble = satZList -> get(j)
        latJavaFloat = satLatList -> get(j)
        lonJavaFloat = satLonList -> get(j)
        ltJavaDouble = satLTList -> get(j)

        x[j] = kmToRe(xJavaDouble -> doubleValue())
        y[j] = kmToRe(yJavaDouble -> doubleValue())
        z[j] = kmToRe(zJavaDouble -> doubleValue())
        lat[j] = latJavaFloat -> floatValue()
        lon[j] = lonJavaFloat -> floatValue()
        localTime[j] = ltJavaDouble -> doubleValue()

       satTimeJ = satTime -> get(j)
       ;
       ; note that satTimeJ is of type javax.xml.datatype.XMLGregorianCalendar
       ; refer to the XMLGregorianCalendar documentation for access to
       ; individual date/time values (e.g., year, hour, etc.)
       ;
       print, satTimeJ -> toString(), x[j], y[j], z[j], lat[j], lon[j], $
              localTime[j]

    ENDFOR

;    IPLOT, x, y, z, /FIT_TO_VIEW, OVERPLOT=1, TITLE='Orbit Plot'

ENDFOR 

OBJ_DESTROY, system, satelliteSituationCenterService, sscWsUrl, sscWsQName, ssc
OBJ_DESTROY, datatypeFactory, dataRequest


END
