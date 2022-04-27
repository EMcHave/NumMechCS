
namespace NumMechCS
{
    class Program
    {
        static void Main(string[] args)
        {
            FEMThermal heattask = new FEMThermal("heatdammesh.txt", "heatdamBC.txt",
                new Concrete());
            heattask.Solve();
            heattask.writeVKT();
            heattask.tempComparation();
        }
    }
}