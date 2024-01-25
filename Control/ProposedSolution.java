// STAR-CCM+ macro: ProposedSolution.java
// Written by STAR-CCM+ 11.06.011
package macro;

import java.util.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import star.common.*;
import star.base.neo.*;
import star.base.report.*;
import star.flow.*;

public class ProposedSolution extends StarMacro {

    public void execute() {
        execute0();
    }

    private void execute0() {

        Simulation simulation_0 = getActiveSimulation();

        Solution solution_0 = simulation_0.getSolution();

        solution_0.clearSolution(); // clears any data from past simulations
        // getting the report for each of the probes
        solution_0.initializeSolution();

        // Everything here is to import matrices into star ccm

        String csvFile = "/home/wh219/FYP/Controller/A2_9000_2.csv";
        String line = "";
        String csvSplitBy = ",";

        double[][] matrixA = new double[50][50]; // initialize matrix with appropriate dimensions

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {

            int i = 0; // initialize row counter

            while ((line = br.readLine()) != null) {

                String[] row = line.split(csvSplitBy); // split line into array of values

                for (int j = 0; j < row.length; j++) {
                    matrixA[i][j] = Double.parseDouble(row[j]); // convert value to double and store in matrix
                }

                i++; // increment row counter

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        String csvFile2 = "/home/wh219/FYP/Controller/B2_9000_2.csv";

        double[][] matrixB = new double[50][4]; // initialize matrix with appropriate dimensions

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile2))) {

            int i = 0; // initialize row counter

            while ((line = br.readLine()) != null) {

                String[] row = line.split(csvSplitBy); // split line into array of values

                for (int j = 0; j < row.length; j++) {
                    matrixB[i][j] = Double.parseDouble(row[j]); // convert value to double and store in matrix
                }

                i++; // increment row counter

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        String csvFile3 = "/home/wh219/FYP/Controller/C2_9000_2.csv";

        double[][] matrixC = new double[50][50]; // initialize matrix with appropriate dimensions

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile3))) {

            int i = 0; // initialize row counter

            while ((line = br.readLine()) != null) {

                String[] row = line.split(csvSplitBy); // split line into array of values

                for (int j = 0; j < row.length; j++) {
                    matrixC[i][j] = Double.parseDouble(row[j]); // convert value to double and store in matrix
                }

                i++; // increment row counter

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        String csvFile4 = "/home/wh219/FYP/Controller/G2_9000_2.csv";

        double[][] matrixG = new double[20][50]; // initialize matrix with appropriate dimensions

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile4))) {

            int i = 0; // initialize row counter

            while ((line = br.readLine()) != null) {

                String[] row = line.split(csvSplitBy); // split line into array of values

                for (int j = 0; j < row.length; j++) {
                    matrixG[i][j] = Double.parseDouble(row[j]); // convert value to double and store in matrix
                }

                i++; // increment row counter

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        String csvFile5 = "/home/wh219/FYP/Controller/K2_9000_2.csv";

        double[][] matrixK = new double[4][50]; // initialize matrix with appropriate dimensions

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile5))) {

            int i = 0; // initialize row counter

            while ((line = br.readLine()) != null) {

                String[] row = line.split(csvSplitBy); // split line into array of values

                for (int j = 0; j < row.length; j++) {
                    matrixK[i][j] = Double.parseDouble(row[j]); // convert value to double and store in matrix
                }

                i++; // increment row counter

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        String csvFile6 = "/home/wh219/FYP/Controller/L2_9000_2.csv";

        double[][] matrixL = new double[50][20]; // initialize matrix with appropriate dimensions

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile6))) {

            int i = 0; // initialize row counter

            while ((line = br.readLine()) != null) {

                String[] row = line.split(csvSplitBy); // split line into array of values

                for (int j = 0; j < row.length; j++) {
                    matrixL[i][j] = Double.parseDouble(row[j]); // convert value to double and store in matrix
                }

                i++; // increment row counter

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        String csvFile7 = "/home/wh219/FYP/Controller/x02_9000_2.csv";

        double[] vectorx0 = new double[50]; // initialize matrix with appropriate dimensions

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile7))) {

            int i = 0; // initialize row counter

            while ((line = br.readLine()) != null) {

                String[] row = line.split(csvSplitBy); // split line into array of values

                for (int j = 0; j < row.length; j++) {
                    vectorx0[i] = Double.parseDouble(row[j]); // convert value to double and store in matrix
                }

                i++; // increment row counter

            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        // put into apache format of a matrix//

        RealMatrix A = new Array2DRowRealMatrix(matrixA);
        RealMatrix B = new Array2DRowRealMatrix(matrixB);
        RealMatrix C = new Array2DRowRealMatrix(matrixC);
        RealMatrix G = new Array2DRowRealMatrix(matrixG);
        RealMatrix K = new Array2DRowRealMatrix(matrixK);
        RealMatrix L = new Array2DRowRealMatrix(matrixL);
        RealVector x0 = new ArrayRealVector(vectorx0);
        RealVector u = new ArrayRealVector(new double[] { 0, 0, 0, 0 }, false); // initialise no actuation

        // Solution solution_0 = simulation_0.getSolution();

        // solution_0.clearSolution(); // clears any data from past simulations
        // getting the report for each of the probes
        // solution_0.initializeSolution();

        MaxReport maxReport_1 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 1"));

        MaxReport maxReport_2 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 2"));

        MaxReport maxReport_3 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 3"));

        MaxReport maxReport_4 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 4"));

        MaxReport maxReport_5 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 5"));

        MaxReport maxReport_6 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 6"));

        MaxReport maxReport_0 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 7"));

        MaxReport maxReport_7 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 8"));

        MaxReport maxReport_8 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 9"));

        MaxReport maxReport_9 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 10"));

        MaxReport maxReport_11 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 11"));

        MaxReport maxReport_12 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 12"));

        MaxReport maxReport_13 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 13"));

        MaxReport maxReport_14 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 14"));

        MaxReport maxReport_15 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 15"));

        MaxReport maxReport_16 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 16"));

        MaxReport maxReport_17 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 17"));

        MaxReport maxReport_18 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 18"));

        MaxReport maxReport_19 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 19"));

        MaxReport maxReport_20 = ((MaxReport) simulation_0.getReportManager().getReport("Probe 20"));

        //

        Region region_0 = simulation_0.getRegionManager().getRegion("Body 1");

        // Boundary boundary_1 =
        // region_0.getBoundaryManager().getBoundary("bottom_top_cylinder");
        Boundary boundary_0 = region_0.getBoundaryManager().getBoundary("Actuator 1");
        Boundary boundary_1 = region_0.getBoundaryManager().getBoundary("Actuator 2");
        Boundary boundary_2 = region_0.getBoundaryManager().getBoundary("Actuator 3");
        Boundary boundary_3 = region_0.getBoundaryManager().getBoundary("Actuator 4");

        VelocityMagnitudeProfile velocityMagnitudeProfile_0 = boundary_0.getValues()
                .get(VelocityMagnitudeProfile.class);

        velocityMagnitudeProfile_0.setMethod(ConstantScalarProfileMethod.class);

        VelocityMagnitudeProfile velocityMagnitudeProfile_1 = boundary_1.getValues()
                .get(VelocityMagnitudeProfile.class);

        velocityMagnitudeProfile_1.setMethod(ConstantScalarProfileMethod.class);

        VelocityMagnitudeProfile velocityMagnitudeProfile_2 = boundary_2.getValues()
                .get(VelocityMagnitudeProfile.class);

        velocityMagnitudeProfile_2.setMethod(ConstantScalarProfileMethod.class);

        VelocityMagnitudeProfile velocityMagnitudeProfile_3 = boundary_3.getValues()
                .get(VelocityMagnitudeProfile.class);

        velocityMagnitudeProfile_3.setMethod(ConstantScalarProfileMethod.class);

        velocityMagnitudeProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(0);
        velocityMagnitudeProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(0);
        velocityMagnitudeProfile_2.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(0);
        velocityMagnitudeProfile_3.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(0);

        simulation_0.getSimulationIterator().step(10000); // do 10000 timesteps before implementing controller 10000 =
                                                          // 100 seconds

        // Loop over 1000 iterations
        for (int i = 0; i < 9000; i++) {

            // Take measurements of velocity from STAR-CCM+
            double s1 = maxReport_1.getReportMonitorValue();
            double s2 = maxReport_2.getReportMonitorValue();
            double s3 = maxReport_3.getReportMonitorValue();
            double s4 = maxReport_4.getReportMonitorValue();
            double s5 = maxReport_5.getReportMonitorValue();
            double s6 = maxReport_6.getReportMonitorValue();
            double s7 = maxReport_0.getReportMonitorValue(); // this is u velocity
            double s8 = maxReport_7.getReportMonitorValue();
            double s9 = maxReport_8.getReportMonitorValue();
            double s10 = maxReport_9.getReportMonitorValue();
            double s11 = maxReport_11.getReportMonitorValue();
            double s12 = maxReport_12.getReportMonitorValue();
            double s13 = maxReport_13.getReportMonitorValue();
            double s14 = maxReport_14.getReportMonitorValue();
            double s15 = maxReport_15.getReportMonitorValue();
            double s16 = maxReport_16.getReportMonitorValue();
            double s17 = maxReport_17.getReportMonitorValue(); // this is u velocity
            double s18 = maxReport_18.getReportMonitorValue();
            double s19 = maxReport_19.getReportMonitorValue();
            double s20 = maxReport_20.getReportMonitorValue();

            RealVector S = new ArrayRealVector(new double[] { s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13,
                    s14, s15, s16, s17, s18, s19, s20 }, false);

            ///////////////////////////////////////////////

            // x0 = A * x0 + L * (s - G * x0) + B * u; // this is my estimator equation
            RealVector x_e = A.operate(x0).add(L.operate(S.subtract(G.operate(x0)))).add(B.operate(u));

            // now calculate the control law
            u = K.operate(x_e.mapMultiply(-1));
            // u = -K * x0;

            double u1 = u.getEntry(0);
            double u2 = u.getEntry(1);
            double u3 = u.getEntry(2);
            double u4 = u.getEntry(3);

            // assign the actuation into the boundaries

            velocityMagnitudeProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(u1);
            velocityMagnitudeProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(u2);
            velocityMagnitudeProfile_2.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(u3);
            velocityMagnitudeProfile_3.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(u4);

            // set x_e to be x0 for the next loop
            x0 = x_e;
            // Run 1 Step
            // simulation_0.println("A matrix = " + matrixA[1][2]); // 0.01798
            simulation_0.println("K matrix = " + K.getEntry(1, 2)); // 0.009155
            simulation_0.println("1st actuation = " + u1);
            simulation_0.println("2nd actuation = " + u2);
            simulation_0.println("3rd actuation = " + u3);
            simulation_0.println("4th actuation = " + u4);
            // System.out.println("x0: v2 " + vectorx0[i]);

            // this should match the extraction rate so if timestep is 0.001 but data
            // extraction rate is at 0.01 then controller should be the same as 0.01 so have
            // it step by 10 i.e 10*0.001 = 0.01
            simulation_0.getSimulationIterator().step(10);

        }
        // Print last achieved lift
        // simulation_0.println("Lift achieved = " +
        // forceReport_0.getReportMonitorValue());

        // DONE
    }
}