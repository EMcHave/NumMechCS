using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumMechCS
{
    internal struct Node //описывает узел
    {
        public int id; //номер узла
        public double x;
        public double y; //координаты
    }

    internal struct Constraint //закрепление
    {
        public int nodeId; //номер закрепленного узла
        public bool isXfixed; //закреплено ли смещение по x
        public bool isYfixed; //по y
        public bool isAngleFixed; //закреплен ли поворот
    }

    internal struct CForce //сосредоточенная сила и ее координаты
    {
        public int nodeId; // к какому узлу приложена
        public float Fx;
        public float Fy;
    }

    internal struct SForce //распределенная нагрузка
    {
        public List<Node> nodes; //массив узлов, к которым она приложена приложена
        public double[] StartEndMultiplier;
        public float Fx;
        public float Fy;
        public double loadLength
        {
            get
            {
                double l2 = Math.Pow(nodes[nodes.Count].x - nodes[0].x, 2) +
                     Math.Pow(nodes[nodes.Count].y - nodes[0].y, 2);
                return Math.Sqrt(l2);
            }
        }
        public double loadMultiplier
        {
            get
            {
                return (StartEndMultiplier[1] - StartEndMultiplier[0]) / loadLength; ;
            }
        }
        public double faceLength(int i)
        {
            double l2 = Math.Pow(nodes[i+1].x - nodes[i].x, 2) +
                     Math.Pow(nodes[i+1].y - nodes[i].y, 2);
            return Math.Sqrt(l2);
        }
    }

    interface IMaterial
    {
        double E { get; } //модуль Юнга
        double V { get; } //коэф пуассона
        double ro { get; } //плотность в кг/м^3
    }
}
