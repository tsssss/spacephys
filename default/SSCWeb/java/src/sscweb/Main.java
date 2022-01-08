/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package sscweb;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Date;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;
import java.util.SimpleTimeZone;
import java.text.SimpleDateFormat;

import javax.xml.datatype.DatatypeFactory;
import javax.xml.datatype.XMLGregorianCalendar;
import javax.xml.namespace.QName;
import javax.xml.ws.BindingProvider;

import com.sun.xml.ws.developer.JAXWSProperties;
import gov.nasa.gsfc.spdf.ssc.client.*;

/**
 *
 * @author Sheng
 */
public class Main {

    private static final SimpleTimeZone UTC_TIME_ZONE = new SimpleTimeZone(0,"UTC");
    private static final SimpleDateFormat DATE_FORMATTER = new SimpleDateFormat("yyyy/DDD HH:mm:ss");
    static {DATE_FORMATTER.setTimeZone(UTC_TIME_ZONE);}
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args)
    throws Exception {

        System.setProperty("http.agent", "WsExample (" + System.getProperty("os.name") + " " + System.getProperty("os.arch") + ")");
        String WSDL_URL = "http://sscweb.gsfc.nasa.gov/WS/ssc/2/SatelliteSituationCenterService?wsdl";
        SatelliteSituationCenterService service = new SatelliteSituationCenterService(
            new URL(args[0]),
            new QName("http://ssc.spdf.gsfc.nasa.gov/", "SatelliteSituationCenterService"));
        SatelliteSituationCenterInterface ssc = service.getSatelliteSituationCenterPort();
        Map<String, Object> requestContextMap = ((BindingProvider)ssc).getRequestContext();
        String dumpMsgs = System.getProperty("com.sun.xml.ws.transport.http.client.HttpTransportPipe.dump");

    }


}
