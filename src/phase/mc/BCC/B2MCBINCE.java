/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package phase.mc.BCC;

import io.prnt;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author metallurgy
 */
public class B2MCBINCE extends BCCMCBINCE {
    //phase specific parameters
    private String phaseTag_local = "B2MCBINCE";
    private int lc_local[] = {1, 2, 2, 1, 1};
    private int tc_local = 7;
    private int nxc_local = 1;
    private int nc_local = 5;

    private int lcf_local[] = {1, 2, 2, 1, 1};
    private int tcf_local = 7;
    private int nxcf_local = 1;
    private int ncf_local = 5;

    private int rcf_local[] = {2, 2, 3, 4, 1};

    private double m_local[][] = {
        {1},
        {0.5, 0.5}, // αα and ββ second-neighbour correlations
        {0.5, 0.5}, // ααβ and αββ triangle correlations
        {1},        // tetrahedron
        {1}         // point
    };

    private String mcMethod_local = "exchange";

    public B2MCBINCE(double[] eci, double T, double xB, int latticeSize, int latticeType) throws IOException {
        //phase specific parameters
        prnt.writeln("  B2MCBINCE constructor method called");
        setPhaseTag(phaseTag_local);
        setLc(lc_local);
        setTc(tc_local);
        setNxc(nxc_local);
        setNc(nc_local);
        setLcf(lcf_local);
        setTcf(tcf_local);
        setNxcf(nxcf_local);
        setNcf(ncf_local);
        setM(m_local);
        //macroscopic parameters
        setR(((Math.abs(eci[0]) == 1) ? 1 : 8.3144)); //  Universal gas constant
        setEdis(eci);
        setT(T);
        setX(xB);
        //MC specfic parameters
        setMcMethod(mcMethod_local);
        setlatticeSize(latticeSize);
        setNLP(calNLP(latticeSize));
        setConfig(genConfig(xB, latticeSize, latticeType));

        double eta = calLongRangeOrder();

        if(Math.abs(xB - 0.5) < 1e-8 && eta < 0.99){
            throw new RuntimeException(
                "Perfect B2 initialization failed. eta = " + eta);
        }

        setInitU(); // calculate correlation functions

        // System.out.println("--------------------------------");
        // prnt.writeln("InitialLRO:" + calLongRangeOrder());
        // System.out.println("--------------------------------");

        prnt.writeln("B2MCBINCE constructor method ended");
    }
//Getter Methods
//Setter  Methods
//public methods
//private methods

    private int[][][] genConfig(double X, int latticeSizeIn, int latticeType) throws IOException {
        prnt.writeln("B2MCBINCE.genConfig() called");
        int nlp = calNLP(latticeSizeIn);
        int[][][] iconfig = new int[latticeSizeIn][latticeSizeIn][latticeSizeIn];
        int[][] coordArray = coordstore(latticeSizeIn);
        ArrayList<int[]> cornerSites = new ArrayList<>();
        ArrayList<int[]> bodySites   = new ArrayList<>();

       
        // new addtion
        int totalB = (int)(X * nlp);

    
    
        int BonBodyCenter = 0;
        int BonCorner = 0;

        // int as = 0;
        // int bs = 0;
        // int no_B = (int) (X * nlp);//no of B element.
        // if (latticeType == 1) {
        //     //randomly call bc and corner sites and assign them A or B atoms
        //     for (int l1 = 1; l1 <= nlp; l1++) {
        //         int rand = (int) (Math.random() * (nlp - l1));
        //         if (l1 <= no_B) {
        //             iconfig[coordArray[rand][0]][coordArray[rand][1]][coordArray[rand][2]] = +1;
        //             bs++;
        //         } else {
        //             iconfig[coordArray[rand][0]][coordArray[rand][1]][coordArray[rand][2]] = -1;
        //             as++;
        //         }
        //         coordArray[rand][0] = coordArray[nlp - l1][0];
        //         coordArray[rand][1] = coordArray[nlp - l1][1];
        //         coordArray[rand][2] = coordArray[nlp - l1][2];
        //     }
        //     System.out.println("X=" + X + ", filled B Atoms=" + bs);
        // }
       

            for (int i = 0; i < coordArray.length; i++) {

                int x = coordArray[i][0];
                int y = coordArray[i][1];
                int z = coordArray[i][2];
            
                if ((x % 2 == 0) &&
                    (y % 2 == 0) &&
                    (z % 2 == 0)) {
            
                    iconfig[x][y][z] = -1;
            
                    cornerSites.add(coordArray[i]);
            
                } else {
            
                    iconfig[x][y][z] = 1;
            
                    bodySites.add(coordArray[i]);
            
                    BonBodyCenter++;
                }
            }
        final int currentB = BonBodyCenter;
        Collections.shuffle(cornerSites);
        Collections.shuffle(bodySites);

        // ---------- Case 1 : xB > 0.5 ----------

        if (totalB > currentB) {

            int extraB = totalB - currentB;

            if (extraB > cornerSites.size()) {
                throw new IllegalArgumentException("Composition exceeds available corner sites.");
            }

            for (int i = 0; i < extraB; i++) {

                int[] c = cornerSites.get(i);
            
                iconfig[c[0]][c[1]][c[2]] = 1;
            
                BonCorner++;
            }
        }
        // ---------- Case 2 : xB < 0.5 ----------

        if (totalB < currentB) {

            int removeB = currentB - totalB;

            if (removeB > bodySites.size()) {
                throw new IllegalArgumentException("Too few body-center B atoms.");
            }

            for (int i = 0; i < removeB; i++) {

                int[] c = bodySites.get(i);
            
                iconfig[c[0]][c[1]][c[2]] = -1;
            
                BonBodyCenter--;
            }
        }
        System.out.println("--------------------------------");

        System.out.println("Requested xB        = " + X);
        
        System.out.println("Total lattice sites = " + nlp);
        
        System.out.println("Required B atoms    = " + totalB);
        
        System.out.println("Body-center B       = " + BonBodyCenter);
        
        System.out.println("Corner B            = " + BonCorner);
        
        System.out.println("Actual B atoms      = " + (BonBodyCenter + BonCorner));
        
        System.out.println("--------------------------------");

        prnt.writeln("B2MCBINCE.genConfig() ended");
        return (iconfig);
    }

    private void setInitU() throws IOException {//Phase specific method
        double[] u_In = calU();
        setU(u_In);
    }

    @Override
    public double[] calU() throws IOException {
        prnt.writeln("B2MCBINCE.calU() method called", 1);
        int[][][] config = getConfig();
        int size = getLatticeSize();
        int nlp = calNLP(size);
        int tcf = getTcf();
        int siteOperators[] = new int[tcf];
        double[] corrFunc = new double[tcf];
        for (int temp = 0; temp < tcf; temp++) {
            siteOperators[temp] = 0;
        }
        for (int i = 0; i < size; i = i + 1) {
            for (int j = 0; j < size; j = j + 1) {
                for (int k = 0; k < size; k = k + 1) {
                    if (((i % 2 == 0) && (j % 2 == 0) && (k % 2 == 0)) || ((i % 2 == 1) && (j % 2 == 1) && (k % 2 == 1))) {//Visiting all BCC sites
                        int site = config[i][j][k];
                        int[] nnsites = firstNeighbourSites(config, size, i, j, k);
                        int[] nsites = secondNeighbourSites(config, size, i, j, k);
                        //prnt.list(nnsites, "nnsites");
                        siteOperators[0] = siteOperators[0] + site * (nnsites[0] + nnsites[1] + nnsites[2] + nnsites[3] + nnsites[4] + nnsites[5] + nnsites[6] + nnsites[7]);
                        if ((i % 2 == 0) && (j % 2 == 0) && (k % 2 == 0)) {

                            // alpha sublattice: alpha-alpha 2NN correlation
                            siteOperators[1] = siteOperators[1]
                                    + site * (nsites[0] + nsites[1] + nsites[2]
                                    + nsites[3] + nsites[4] + nsites[5]);
                        
                        } else {
                        
                            // beta sublattice: beta-beta 2NN correlation
                            siteOperators[2] = siteOperators[2]
                                    + site * (nsites[0] + nsites[1] + nsites[2]
                                    + nsites[3] + nsites[4] + nsites[5]);
                        }
                        // Triangle contributions from first-neighbour pairs
                        int triangleFrom1NN =
                        (nnsites[0] * nnsites[1]
                        + nnsites[1] * nnsites[3]
                        + nnsites[3] * nnsites[2]
                        + nnsites[2] * nnsites[0])
                      + (nnsites[4] * nnsites[5]
                        + nnsites[5] * nnsites[7]
                        + nnsites[7] * nnsites[6]
                        + nnsites[6] * nnsites[4])
                      + (nnsites[1] * nnsites[5]
                        + nnsites[4] * nnsites[0])
                      + (nnsites[3] * nnsites[7]
                        + nnsites[6] * nnsites[2]);
                
                // Triangle contributions involving second neighbours
                        int triangleFrom2NN =
                                nsites[0] * (nnsites[0] + nnsites[1]
                                            + nnsites[3] + nnsites[2])
                            + nsites[1] * (nnsites[4] + nnsites[5]
                                            + nnsites[7] + nnsites[6])
                            + nsites[2] * (nnsites[0] + nnsites[1]
                                            + nnsites[5] + nnsites[4])
                            + nsites[3] * (nnsites[2] + nnsites[3]
                                            + nnsites[7] + nnsites[6])
                            + nsites[4] * (nnsites[0] + nnsites[2]
                                            + nnsites[6] + nnsites[4])
                            + nsites[5] * (nnsites[1] + nnsites[3]
                                            + nnsites[7] + nnsites[5]);
                        
                        // Alpha sublattice: even-even-even coordinates
                        if ((i % 2 == 0) && (j % 2 == 0) && (k % 2 == 0)) {
                        
                            siteOperators[3] += site * triangleFrom2NN; // ααβ
                            siteOperators[4] += site * triangleFrom1NN; // αββ
                        
                        } else {
                        
                            siteOperators[3] += site * triangleFrom1NN; // ααβ
                            siteOperators[4] += site * triangleFrom2NN; // αββ
                        }
                        siteOperators[5] = siteOperators[5] + site * (nsites[0] * (nnsites[0] * nnsites[1] + nnsites[1] * nnsites[3] + nnsites[3] * nnsites[2] + nnsites[2] * nnsites[0])
                                + nsites[1] * (nnsites[4] * nnsites[5] + nnsites[5] * nnsites[7] + nnsites[7] * nnsites[6] + nnsites[6] * nnsites[4])
                                + nsites[2] * (nnsites[0] * nnsites[1] + nnsites[1] * nnsites[5] + nnsites[5] * nnsites[4] + nnsites[4] * nnsites[0])
                                + nsites[3] * (nnsites[2] * nnsites[3] + nnsites[3] * nnsites[7] + nnsites[7] * nnsites[6] + nnsites[6] * nnsites[2])
                                + nsites[4] * (nnsites[0] * nnsites[2] + nnsites[2] * nnsites[6] + nnsites[6] * nnsites[4] + nnsites[4] * nnsites[0])
                                + nsites[5] * (nnsites[1] * nnsites[3] + nnsites[3] * nnsites[7] + nnsites[7] * nnsites[5] + nnsites[5] * nnsites[1]));
                        siteOperators[6] = siteOperators[6] + site;
                    }
                }
            }
        }
        int[] parentCluster = {0, 1, 1, 2, 2, 3, 4};

        for (int i = 0; i < tcf; i++) {
            double denominator =
                    getRcdis()[parentCluster[i]]
                    * getMgdis()[parentCluster[i]]
                    * nlp;
                    if (i == 1 || i == 2 || i == 3 || i == 4) {
                        denominator *= 0.5;
                    }
            corrFunc[i] =
                    (1.0 * siteOperators[i]) / denominator;
            // System.out.println("siteOperators[" + i + "]=" + siteOperators[i]);
            // System.out.println("corrFunc[" + i + "]=" + corrFunc[i]);
        }
        prnt.writeln("B2MCBINCE.calU() method ended", 1);
        return (corrFunc);
    }
    public double calLongRangeOrder() {

        int[][][] config = getConfig();
        int size = getLatticeSize();
    
        int cornerA = 0;
        int bodyB = 0;
    
        int totalCorner = 0;
        int totalBody = 0;
    
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                for (int k = 0; k < size; k++) {
    
                    // Visit only valid BCC sites
                    if (((i % 2 == 0) && (j % 2 == 0) && (k % 2 == 0)) ||
                        ((i % 2 == 1) && (j % 2 == 1) && (k % 2 == 1))) {
    
                        int site = config[i][j][k];
    
                        if ((i % 2 == 0) &&
                            (j % 2 == 0) &&
                            (k % 2 == 0)) {
    
                            totalCorner++;
    
                            if (site == -1)
                                cornerA++;
    
                        } else {
    
                            totalBody++;
    
                            if (site == 1)
                                bodyB++;
                        }
                    }
                }
            }
        }
    
        double etaCorner = (double) cornerA / totalCorner;
        double etaBody   = (double) bodyB / totalBody;
    
        double eta = etaCorner + etaBody - 1.0;

        // Numerical safety
        if (eta > 1.0) eta = 1.0;
        if (eta < 0.0) eta = 0.0;

        return eta;
            }


            @Override
public double calLHc(int[] siteCoord) throws IOException {

    int[][][] config = getConfig();
    int size = getLatticeSize();
    double[] edis = getEcdis();

    int[] SiteCF = new int[7];

    int site = config[siteCoord[0]][siteCoord[1]][siteCoord[2]];

    int[] nnsites = firstNeighbourSites(
            config, size,
            siteCoord[0], siteCoord[1], siteCoord[2]);

    int[] nsites = secondNeighbourSites(
            config, size,
            siteCoord[0], siteCoord[1], siteCoord[2]);

    // CF[0] : 1NN alpha-beta
    SiteCF[0] = site * (
            nnsites[0] + nnsites[1] + nnsites[2] + nnsites[3]
            + nnsites[4] + nnsites[5] + nnsites[6] + nnsites[7]
    );

    // 2NN contribution
    int secondNN = site * (
            nsites[0] + nsites[1] + nsites[2]
            + nsites[3] + nsites[4] + nsites[5]
    );

    // CF[1] : alpha-alpha 2NN
    // CF[2] : beta-beta 2NN
    if ((siteCoord[0] % 2 == 0)
            && (siteCoord[1] % 2 == 0)
            && (siteCoord[2] % 2 == 0)) {

        SiteCF[1] = secondNN;

    } else {

        SiteCF[2] = secondNN;
    }

    // CF[3] : triplet
        int triangleFrom1NN =
            (nnsites[0] * nnsites[1]
            + nnsites[1] * nnsites[3]
            + nnsites[3] * nnsites[2]
            + nnsites[2] * nnsites[0])
            + (nnsites[4] * nnsites[5]
            + nnsites[5] * nnsites[7]
            + nnsites[7] * nnsites[6]
            + nnsites[6] * nnsites[4])
            + (nnsites[1] * nnsites[5]
            + nnsites[4] * nnsites[0])
            + (nnsites[3] * nnsites[7]
            + nnsites[6] * nnsites[2]);

        int triangleFrom2NN =
            nsites[0] * (nnsites[0] + nnsites[1]
                        + nnsites[3] + nnsites[2])
            + nsites[1] * (nnsites[4] + nnsites[5]
                        + nnsites[7] + nnsites[6])
            + nsites[2] * (nnsites[0] + nnsites[1]
                        + nnsites[5] + nnsites[4])
            + nsites[3] * (nnsites[2] + nnsites[3]
                        + nnsites[7] + nnsites[6])
            + nsites[4] * (nnsites[0] + nnsites[2]
                        + nnsites[6] + nnsites[4])
            + nsites[5] * (nnsites[1] + nnsites[3]
                        + nnsites[7] + nnsites[5]);

        // CF[3] = ααβ; CF[4] = αββ
        if ((siteCoord[0] % 2 == 0)
            && (siteCoord[1] % 2 == 0)
            && (siteCoord[2] % 2 == 0)) {

        SiteCF[3] = site * triangleFrom2NN;
        SiteCF[4] = site * triangleFrom1NN;

        } else {

        SiteCF[3] = site * triangleFrom1NN;
        SiteCF[4] = site * triangleFrom2NN;
        }

    // CF[5] : quadruplet
    SiteCF[5] =
            site * (
                nsites[0] * (
                    nnsites[0] * nnsites[1]
                    + nnsites[1] * nnsites[3]
                    + nnsites[3] * nnsites[2]
                    + nnsites[2] * nnsites[0])

                + nsites[1] * (
                    nnsites[4] * nnsites[5]
                    + nnsites[5] * nnsites[7]
                    + nnsites[7] * nnsites[6]
                    + nnsites[6] * nnsites[4])

                + nsites[2] * (
                    nnsites[0] * nnsites[1]
                    + nnsites[1] * nnsites[5]
                    + nnsites[5] * nnsites[4]
                    + nnsites[4] * nnsites[0])

                + nsites[3] * (
                    nnsites[2] * nnsites[3]
                    + nnsites[3] * nnsites[7]
                    + nnsites[7] * nnsites[6]
                    + nnsites[6] * nnsites[2])

                + nsites[4] * (
                    nnsites[0] * nnsites[2]
                    + nnsites[2] * nnsites[6]
                    + nnsites[6] * nnsites[4]
                    + nnsites[4] * nnsites[0])

                + nsites[5] * (
                    nnsites[1] * nnsites[3]
                    + nnsites[3] * nnsites[7]
                    + nnsites[7] * nnsites[5]
                    + nnsites[5] * nnsites[1])
            );

    // CF[5] : point correlation
    SiteCF[6] = site;

    /*
     * Five parent clusters:
     *
     * ECI[0] -> CF[0]
     * ECI[1] -> 0.5*CF[1] + 0.5*CF[2]
     * ECI[2] -> CF[3]
     * ECI[3] -> CF[4]
     * ECI[4] -> CF[5]
     *
     * This matches:
     * m_local = {{1},{1},{0.5,0.5},{1},{1}}
     */
    double SiteEnthalpy =
        edis[0] * SiteCF[0]
      + edis[1] * (0.5 * SiteCF[1] + 0.5 * SiteCF[2])
      + edis[2] * (0.5 * SiteCF[3] + 0.5 * SiteCF[4])
      + edis[3] * SiteCF[5]
      + edis[4] * SiteCF[6];

return SiteEnthalpy;
        }
}
