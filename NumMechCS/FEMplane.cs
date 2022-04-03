using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace NumMechCS
{
    internal class FEMplane 
    {
        
        private Matrix<double>? StiffnessGlobal;
        private Vector<double>? ForcesVector;
        private List<Constraint>? constraints; 
        private List<CForce> concentratedForces; 
        private List<SForce> surfaceFoces;
        private IMaterial material;
        const double g = -9.81;

        public static Matrix<double> Dmatrix;
        
        public List<Node>? nodes { get; } 
        public List<Element>? elements { get; }
        public FEMplane(string inputFileAdress, IMaterial material,
                        bool cf, bool sf, bool gf)
        {
            this.material = material;
            readInputFile(inputFileAdress, material, cf, sf, gf);
            StiffnessGlobal = Matrix<double>.Build.Sparse(2 * nodes.Count(), 2 * nodes.Count());
            ForcesVector = Vector<double>.Build.Sparse(2 * nodes.Count());
            Dmatrix = material.E/(1-Math.Pow(material.V,2))*Matrix<double>.Build.DenseOfArray(
                new double[,] 
                { 
                    { 1, material.V, 0},
                    {material.V, 1, 0},
                    {0,0, (1-material.V)/2 }
                });
            buildForcesVector(cf, sf, gf);
        }

        private Vector<double> displacements { get; set; } 
        private Vector<double> strains { get; set; } 
        private Vector<double> stresses { get; set; } 

        public void Solve(bool calculateStrains, bool calculateStresses)
        {
            buildStiffnessMatrix();
            applyConstraints();

            displacements = StiffnessGlobal.Solve(ForcesVector);
        }

        private void applyConstraints()
        {
            List<int> indicesToConstraint = new List<int>();
            foreach(var constraint in constraints)
            {
                if (constraint.isXfixed)
                    indicesToConstraint.Add(2 * constraint.nodeId + 0);
                if (constraint.isYfixed)
                    indicesToConstraint.Add(2 * constraint.nodeId + 1);
            }

            foreach (int DOF in indicesToConstraint)
            {
                ForcesVector[DOF] = 0;
                StiffnessGlobal.ClearRow(DOF);
                StiffnessGlobal.ClearColumn(DOF);
                StiffnessGlobal[DOF, DOF] = 1;
            }

            
        }
        private void readInputFile(string InputFileAdress, IMaterial material,
                                    bool cf, bool sf, bool gf)
        {

        }
        private void buildStiffnessMatrix()
        {
            foreach (var element in elements)
            {
                element.CalculateStiffnessMatrix(ref StiffnessGlobal, ref Dmatrix);
            }
        }

        private void buildForcesVector(bool cf, bool sf, bool gf)
        {
            Vector<double> concForces = Vector<double>.Build.Sparse(2 * nodes.Count());
            Vector<double> surfForces = Vector<double>.Build.Sparse(2 * nodes.Count());
            Vector<double> gForces = Vector<double>.Build.Sparse(2 * nodes.Count());

            if (cf)
            {
                foreach(var force in concentratedForces)
                {
                    concForces[2 * force.nodeId] = force.Fx;
                    concForces[2 * force.nodeId + 1] = force.Fy;
                }
            }
            if(sf)
            {
                foreach (var load in surfaceFoces)
                {
                    double x1 = 0;
                    for (int i =0; i < nodes.Count; i++)
                    {
                        if (i == 0)
                        {
                            double xl1 = x1;
                            double xl2 = x1 + load.faceLength(i)/2;
                            double S = load.loadMultiplier * (xl1 + xl2) / 2*(load.faceLength(i)/2);
                            surfForces[2 * i] += load.Fx * S;
                            surfForces[2 * i + 1] += load.Fy * S;
                            x1 = xl2;
                        }
                        else if(i == nodes.Count - 1)
                        {
                            double xl1 = x1;
                            double xl2 = x1 + load.faceLength(i-1)/2;
                            double S = load.loadMultiplier * (xl1 + xl2) / 2 * (load.faceLength(i-1)/2);
                            surfForces[2 * i] += load.Fx * S;
                            surfForces[2 * i + 1] += load.Fy * S;
                        }
                        else
                        {
                            double xl1 = x1;
                            double xl2 = x1 + (load.faceLength(i-1) + load.faceLength(i)) / 2;
                            double S = load.loadMultiplier * (xl1 + xl2) / 2 * (load.faceLength(i) + load.faceLength(i-1))/2;
                            surfForces[2 * i] += load.Fx * S;
                            surfForces[2 * i + 1] += load.Fy * S;
                            x1 = xl2;
                        }
                    }
                }   
            }
            if(gf)
            {
                foreach(var element in elements)
                {
                    var gVector = element.Jacobian() * material.ro * g / 6 * Vector<double>.Build.DenseOfArray(new double[6] { 0, 1, 0, 1, 0, 1 });
                    gForces[2 * element.nodesIDs[0] + 1] += gVector[1];
                    gForces[2 * element.nodesIDs[1] + 1] += gVector[3];
                    gForces[2 * element.nodesIDs[2] + 1] += gVector[5];
                }
            }

            ForcesVector = concForces + surfForces + gForces;

        }
        /*
        public Vector<double> computeStrains(Vector<double>? stresses)
        {
            //расчет деформаций при необходимости
        }
        public Vector<double> computeStresses(Vector<double>? strains)
        {
            //расчет напряжений при необходимости
        }


        public void writeVTK()
        {
            //генерация VTK-файла
        }
        */
    }
}
