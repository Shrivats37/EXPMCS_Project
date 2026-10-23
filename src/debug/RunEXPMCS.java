/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package debug;

import java.io.IOException;
import io.prnt;
import mcSampler.mcData;
import mcSampler.mcSampler;
import phase.mc.BCC.A2MCBINCE;
import phase.mc.BCC.B2MCBINCE;
import phase.mc.FCC.A1MCBINCE;
import phase.mc.FCC.L10MCBINCE;
import phase.mc.FCC.L12MCBINCE;
import phase.mc.PHASEMCBINCE;

/**
 *
 * @author metallurgy
 */
public class RunEXPMCS {

    public static void main(String[] args) throws IOException {
        long beg = System.currentTimeMillis();
        double xB = 0.5;


        // double T = 1.0; 

        double[] compositions = {
            0.1, 0.2, 0.3, 0.4, 0.45,
            0.5,
            0.55, 0.6, 0.7, 0.8, 0.9
        };


        // for (int i = 0; i < 1; i++) {
        //     debugRunL10(1.73, xB);
        //     xB = xB + 0.01;
        // }

      
       
        double[] temperatures = {
            1.00,1.10,1.20,1.30,1.40,1.50,1.55,1.60,1.65,1.70,1.72,1.74,1.76,1.78,1.80,1.82,1.84,1.86,
            1.88,1.90,1.95,2.00,2.10,2.20,2.40,2.60,3.00,3.50,4.00,5.00,6.00,7,8,9,10,12,15
        };

        for(double T : temperatures)
            {
                // for(int run=1; run<=5; run++)
                // {
                    debugRunB2(T, xB);
                // }
        }

        // for (double xB : compositions) {
        //     System.out.println("\n>>> Starting Simulation for Composition: " + xB);
        
        //     debugRunA2(T, xB);
        // }
        // ----------------------


        long end = System.currentTimeMillis();
        System.out.println("#MC run took " + (double) (end - beg) / 1000 + " sec");
    }

    private static void debugRunA2(double T_In, double x_In) throws IOException {
        int size = 32;
        int latticeType = 1; 
        //double eBCC[] = {-250.0, -166.6666667, 0.0, 25.0, 0.0};
        //double R = 8.3144;
        //double T = 284.0;
        // double eBCC[] = {-1.0, 0.0, 0.0, 0.0, 0.0};
        double eBCC[] = {+1.0, 0.0, 0.0, 0.0, 0.0};
        double R = 1.0;
        double T = T_In;
        double xB = x_In;
        // double T = 6.20;
        // double xB = 0.5;
        A2MCBINCE a2mcbince = new A2MCBINCE(eBCC, T, xB, size, latticeType);


        a2mcbince.printPhaseInfo();
        mcSampler mcsampler = new mcSampler(4000, 1000, a2mcbince);
        mcsampler.printSamplerInfo();
        mcData mcdata = mcsampler.runMC();
        mcdata.calStat();
        mcdata.printStat(true);

         // -------------------------
        // MACHINE READABLE SUMMARY
        // -------------------------

        System.out.println("=====SIMULATION_RESULT_START=====");

        System.out.println("Temperature=" + T);
        System.out.println("Composition=" + xB);

        System.out.println("InitialLRO=0.0");
        System.out.println("FinalLRO=0.0");

        // leave these for now
        // we'll expose them in Step 2

        System.out.println("=====SIMULATION_RESULT_END=====");
    }

    private static void debugRunB2(double T_In, double x_In) throws IOException {

        int size = 32;
    
        int latticeType = 1;
    
        // double eBCC[] = {-1.0, 0.0, 0.0, 0.0, 0.0};
        double eBCC[] = {+1.0,0,0,0,0};
        // // double[] eBCC = {0.0, 1.0, 0.0, 0.0, 0.0};
        // double[] eBCC = {0.0, -1.0, 0.0, 0.0, 0.0};
    
        double R = 1.0;
    
        double T = T_In;
    
        double xB = x_In;
    
        B2MCBINCE b2mcbince =new B2MCBINCE(eBCC,T,xB,size,latticeType);

        double initialLRO = b2mcbince.calLongRangeOrder();
    
        b2mcbince.printPhaseInfo();
    
        mcSampler mcsampler = new mcSampler(4000, 1000, b2mcbince);
    
        mcsampler.printSamplerInfo();
    
        mcData mcdata = mcsampler.runMC();

        mcdata.calStat();

        mcdata.printStat(true);

        // -------------------------
        // MACHINE READABLE SUMMARY
        // -------------------------

        System.out.println("=====SIMULATION_RESULT_START=====");

        System.out.println("Temperature=" + T);
        System.out.println("Composition=" + xB);

        System.out.println("InitialLRO=" + initialLRO);
        System.out.println("FinalLRO=" + b2mcbince.calLongRangeOrder());

        // leave these for now
        // we'll expose them in Step 2

        System.out.println("=====SIMULATION_RESULT_END=====");
    }

    private static void debugRunA1(double T_In, double x_In) throws IOException {
        int size = 32;
        int latticeType = 1;
        double eFCC[] = {-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double R = 1.0;
        double T = T_In;
        double xB = x_In;
        PHASEMCBINCE a1mcbince = new A1MCBINCE(eFCC, T, xB, size, latticeType);
        a1mcbince.printPhaseInfo();
        mcSampler mcsampler = new mcSampler(4000, 1000, a1mcbince);
        mcsampler.printSamplerInfo();
        mcData mcdata = mcsampler.runMC();
        mcdata.calStat();
        mcdata.printStat(true);
    }

    private static void debugRunL10(double T_In, double x_In) throws IOException {
        int size = 32;
        int latticeType = 2;
        double eFCC[] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double R = 1.0;
        double T = T_In;
        double xB = x_In;
        PHASEMCBINCE l10mcbince = new L10MCBINCE(eFCC, T, xB, size, latticeType);
        l10mcbince.printPhaseInfo();
        mcSampler mcsampler = new mcSampler(4000, 1000, l10mcbince);
        mcsampler.printSamplerInfo();
        mcData mcdata = mcsampler.runMC();
        mcdata.calStat();
        mcdata.printStat(true);
    }

    private static void debugRunL12(double T_In, double x_In) throws IOException {
        int size = 32;
        int latticeType = 1;
        double eFCC[] = {1.0, -0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double R = 1.0;
        double T = T_In;
        double xB = x_In;
        PHASEMCBINCE l12mcbince = new L12MCBINCE(eFCC, T, xB, size, latticeType);
        l12mcbince.printPhaseInfo();
        mcSampler mcsampler = new mcSampler(4000, 1000, l12mcbince);
        mcsampler.printSamplerInfo();
        mcData mcdata = mcsampler.runMC();
        mcdata.calStat();
        mcdata.printStat(true);
    }
}
