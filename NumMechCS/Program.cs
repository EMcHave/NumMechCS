
namespace NumMechCS
{
    class Program
    {
        static void Main()
        {

            var Solver = new FEMplane("DamJob.inp",new Concrete(), false, true, true);
            Solver.Solve(false, false);
            Console.Write(Solver.displacements);
        }
    }
}