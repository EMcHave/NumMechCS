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
        
        protected Matrix<double>? StiffnessGlobal;
        protected Vector<double>? ForcesVector;
        private List<Constraint>? constraints; 
        protected List<CForce> concentratedForces; 
        protected List<SForce> surfaceForces;
        protected IMaterial material;
        const double g = -9.81;

        public static Matrix<double> Dmatrix;

        public double thickness;
        public List<Node>? nodes { get; private set; }
        public List<Element>? elements { get; private set; }
        public FEMplane(string inputFileAdress, IMaterial material)
        {
            this.thickness = 1;
            this.material = material;
            readInputFile(inputFileAdress);
            StiffnessGlobal = Matrix<double>.Build.Sparse(2 * nodes.Count(), 2 * nodes.Count());
            ForcesVector = Vector<double>.Build.Sparse(2 * nodes.Count());
        }

        public Vector<double> displacements { get; set; } 
        private Vector<double> strains { get; set; } 
        private Vector<double> stresses { get; set; } 

        public void Solve(bool calculateStrains, bool calculateStresses,
                          bool cf, bool sf, bool gf)
        {
            Dmatrix = material.E / (1 - Math.Pow(material.V, 2)) * Matrix<double>.Build.DenseOfArray(
                new double[,]
                {
                    { 1, material.V, 0},
                    {material.V, 1, 0},
                    {0,0, (1-material.V)/2 }
                });

            buildStiffnessMatrix();
            buildForcesVector(cf, sf, gf);
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
        private void readInputFile(string InputFileAdress)
        {
            nodes = new List<Node>();
            elements = new List<Element>();
            surfaceForces = new List<SForce>();
            concentratedForces = new List<CForce>();
            constraints = new List<Constraint>();
            using(StreamReader sr = new StreamReader(InputFileAdress))
            {
                string? line;
                string[] numbers;
                char[] separators = new char[2] { ' ', ',' };
                line = sr.ReadLine();

                while((line = sr.ReadLine()) != "*Element")
                {
                    numbers = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    Node node = new Node()
                    {
                        id = int.Parse(numbers[0])-1,
                        x = double.Parse(numbers[1], System.Globalization.CultureInfo.InvariantCulture),
                        y = double.Parse(numbers[2], System.Globalization.CultureInfo.InvariantCulture),
                    };
                    nodes.Add(node);
                }
                while((line = sr.ReadLine()) != "*End")
                {
                    numbers = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    elements.Add(new Element(nodes[int.Parse(numbers[1])-1],
                        nodes[int.Parse(numbers[2])-1],
                        nodes[int.Parse(numbers[3])-1],
                        int.Parse(numbers[0])-1));
                }
            }


            ///////////////////////////////////////////////////////
            List<int> loadedNodesNumbers = new List<int>()
            { 4, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 3};
            List<int> constrainedNodesNumbers = new List<int>()
            { 4, 23, 24, 5, 28, 27, 26, 25, 6 };
            List<Node> loadedNodes = new List<Node>();

            foreach (int N in loadedNodesNumbers)
                loadedNodes.Add(nodes[N-1]);


            surfaceForces.Add(new SForce()
            {
                nodes = loadedNodes,
                StartEndMultiplier = new double[2] { 1, 1 },
                Fx = 100000,
                Fy = 0
            });
            concentratedForces.Add(new CForce()
            {
                nodeId = 1,
                Fx = 100000,
                Fy = -30000,
            });

            foreach (int N in constrainedNodesNumbers)
            {
                constraints.Add(new Constraint()
                {
                    nodeId = N - 1,
                    isXfixed = true,
                    isYfixed = true,
                });
            }
            /////////////////////////////////////////////////////

            
        }
        private void buildStiffnessMatrix()
        {
            foreach (var element in elements)
            {
                element.CalculateStiffnessMatrix(ref StiffnessGlobal, ref Dmatrix);
            }
        }

        protected void buildForcesVector(bool cf, bool sf, bool gf)
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
                foreach (var load in surfaceForces)
                {
                    
                    double x1 = 0;
                    double xl1, xl2, S;
                    for (int i = 0; i < load.nodes.Count; i++)
                    {
                        xl1 = x1;
                        if (i == 0)
                        {
                            xl2 = x1 + load.faceLength(i)/2;
                            S = (2*load.StartEndMultiplier[0]+ load.loadMultiplier * (xl1 + xl2)) / 2 * (load.faceLength(i)/2);
                            x1 = xl2;
                        }
                        else if(i == load.nodes.Count - 1)
                        {
                            xl2 = x1 + load.faceLength(i-1)/2;
                            S = (2 * load.StartEndMultiplier[0] + load.loadMultiplier * (xl1 + xl2)) / 2 * (load.faceLength(i-1)/2);
                        }
                        else
                        {
                            xl2 = x1 + (load.faceLength(i-1) + load.faceLength(i)) / 2;
                            S = (2 * load.StartEndMultiplier[0] + load.loadMultiplier * (xl1 + xl2)) / 2 * (xl2 - xl1);
                            x1 = xl2;
                        }
                        surfForces[2 * load.nodes[i].id] += load.Fx * thickness * S;
                        surfForces[2 * load.nodes[i].id + 1] += load.Fy * thickness * S;
                    }
                }   
            }
            if(gf)
            {
                foreach(var element in elements)
                {
                    var gVector = element.Jacobian().Determinant() * material.ro * g * thickness / 6 * Vector<double>.Build.DenseOfArray(new double[6] { 0, 1, 0, 1, 0, 1 });
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
