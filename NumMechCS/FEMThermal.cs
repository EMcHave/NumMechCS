using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using System.Globalization;
using MathNet.Numerics.Data.Text;

namespace NumMechCS
{
    internal class FEMThermal : FEMplane
    {
        private List<BCTemp> BCTemps = new List<BCTemp>();
        public Vector<double> temperatures { get; private set; } 
        public FEMThermal(string inputFileAdress, string BCtemps, IMaterial material)
            : base(inputFileAdress, material)
        {
            StiffnessGlobal = Matrix<double>.Build.Sparse(nodes.Count(), nodes.Count());
            ForcesVector = Vector<double>.Build.Sparse(nodes.Count());
            temperatures = Vector<double>.Build.Sparse(nodes.Count());
            readBC(BCtemps);
        }

        public Vector<double> Solve()
        {

            buildThermalMatrix();
            setBoundaryTemps();
            temperatures = StiffnessGlobal.Solve(ForcesVector);
            return temperatures;
        }

        private void buildThermalMatrix()
        {
            foreach (var element in elements)
            {
                element.CalculateThermalMatrix(ref StiffnessGlobal);
            }
        }
        private void readBC(string BCtemps)
        {
            using (StreamReader sr = new StreamReader(BCtemps))
            {
                string line;
                string[] data;
                char[] separators = new char[2] { ' ', ',' };
                
                while((line = sr.ReadLine()) != "")
                {
                    List<int> temp = new List<int>();
                    data = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    for (int i = 1; i < data.Length; i++)
                        temp.Add(int.Parse(data[i]) - 1);

                    this.BCTemps.Add(new BCTemp
                    {
                        nodes = temp,
                        temp = int.Parse(data[0])
                    });
                }
                while ((line = sr.ReadLine()) != null)
                {
                    data = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    int N1 = int.Parse(data[0])-1;
                    int N2 = int.Parse(data[1])-1;
                    double l = double.Parse(data[2], System.Globalization.CultureInfo.InvariantCulture);
                    for (int i = N1; i < N2; i++)
                        elements[i].lambda = l;
                }
                
            }
        }

        private void setBoundaryTemps()
        {
            foreach(var BC in BCTemps)
            {
                foreach (int n in BC.nodes)
                {
                    StiffnessGlobal.ClearRow(n);
                    StiffnessGlobal[n, n] = 1;
                    ForcesVector[n] = BC.temp;
                }
            }
        }

        public void writeVKT()
        {
            NumberFormatInfo nfi = new NumberFormatInfo();
            nfi.NumberDecimalSeparator = ".";
            string path = "heat.vtk";

            using(StreamWriter sw = new StreamWriter(path, false))
            {
                sw.Write("# vtk DataFile Version 1.0 \n3D triangulation data \nASCII\n\n");
                sw.WriteLine("DATASET POLYDATA");
                sw.WriteLine($"POINTS {nodes.Count} float");
                foreach(var node in nodes)
                    sw.WriteLine($"{node.x.ToString(nfi)}\t {node.y.ToString(nfi)}\t 0");
                sw.WriteLine($"POLYGONS {elements.Count} {4 * elements.Count}");
                foreach (var el in elements)
                    sw.WriteLine($"3\t {el.nodesIDs[0]}\t {el.nodesIDs[1]}\t {el.nodesIDs[2]}");
                sw.WriteLine($"POINT_DATA {nodes.Count}");
                sw.WriteLine("SCALARS Tempreture float 1\nLOOKUP_TABLE default");
                foreach (double t in temperatures)
                    sw.WriteLine($"{t.ToString(nfi)}");
            }
        }

        public void tempComparation()
        {
            string path = "DamTemp.txt";
            using (StreamReader sr = new StreamReader(path))
            {
                string line;
                string[] data;
                char[] separators = new char[2] { ' ', ',' };
                Vector<double> abTemp = Vector<double>.Build.Dense(nodes.Count);
                int k = 0;
                /*
                while ((line = sr.ReadLine()) != null)
                {
                    data = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    abTemp[k++] = (double.Parse(data[1], CultureInfo.InvariantCulture));
                }
                */

                //Vector<double> comparation = temperatures - abTemp;
                //DelimitedWriter.Write("HeatComp.csv", comparation.ToColumnMatrix(), ",");
                DelimitedWriter.Write("HeatC#.csv", temperatures.ToColumnMatrix(), ",");
                //DelimitedWriter.Write("HeatAbaqus.csv", abTemp.ToColumnMatrix(), ",");
            }
        }
    }
}
