using Microsoft.VisualBasic;
using System.Numerics;

class Program()
{
    public const double Epsilon = 1e-6;
    public static void Main(string[] args)
    {
        double x = 0;
        double y = 0;
        double z = 0;

        x = 1;
        y = 1;
        z = 1;

        //x = 0.5;
        //y = 0.8;
        //z = 0.0;

        double reductionFactor = 0.05;
        double r = 10;
        Func<double, double, double, double, double> penaltyFunc = (x, y, z, r) => PenaltyFunction(x, y, z, r);

        var (pointTimeline, functionCalls, cycles) = FastestDescentWithPenalty(x, y, z, r, penaltyFunc, reductionFactor);

        foreach (var point in pointTimeline)
        {
            Console.WriteLine($"x= {point.Item1}, y= {point.Item2}, z= {point.Item3}, lr= {point.Item4}");
        }
        Console.WriteLine($"Function Calls: {functionCalls}, Cycles: {cycles}");
    }

    public static (List<(double, double, double, double)>, int, int) FastestDescentWithPenalty(
        double x, double y, double z, double r, Func<double, double, double, double, double> penaltyFunc, double redFactor=0.1)
    {
        List<(double, double, double, double)> pointTimeLine = new();
        pointTimeLine.Add((x, y, z, 1));
        int cycles = 0;
        int iteration = 0;
        int targetFunctionCalled = 0;

        while (iteration < 500 && r >= 0)
        {
            Console.WriteLine($"Iteration {iteration + 1}, Penalty Multiplier: {r}");
            (double gradX, double gradY, double gradZ) = PartialGradientsWithPenalty(x, y, z, r);
            targetFunctionCalled += 6; // 3 gradinets, each with equality and inequality constraints

            if (Math.Abs(gradX) < Epsilon && Math.Abs(gradY) < Epsilon && Math.Abs(gradZ) < Epsilon)
            {
                Console.WriteLine($"Stopped. gradX: {gradX}, gradY: {gradY}, gradZ: {gradZ}");
                break;
            }

            Func<double, double> goldenSectionFunc = lambda => penaltyFunc(x - lambda * gradX, y - lambda * gradY, z - lambda * gradZ, r);

            double lambda = GoldenSectionSearch(goldenSectionFunc, ref cycles, ref targetFunctionCalled);
            x -= lambda * gradX;
            y -= lambda * gradY;
            z -= lambda * gradZ;

            pointTimeLine.Add((x, y, z, lambda));
            cycles++;
            iteration++;
            r -= redFactor;
        }

        return (pointTimeLine, targetFunctionCalled, cycles);
    }

    private static double GoldenSectionSearch(Func<double, double> func, ref int cycles, ref int targetFunctionCalled, double leftBorder = 0, double rightBorder = 5)
    {
        double T = (-1 + Math.Sqrt(5)) / 2;
        double interval = rightBorder - leftBorder;
        double lambda1 = rightBorder - T * interval;
        double lambda2 = leftBorder + T * interval;

        double f1 = func(lambda1);
        double f2 = func(lambda2);
        targetFunctionCalled += 12;

        while (interval > Epsilon)
        {
            if (f2 < f1)
            {
                leftBorder = lambda1;
                lambda1 = lambda2;
                f1 = f2;
                interval = rightBorder - leftBorder;
                lambda2 = leftBorder + T * interval;
                f2 = func(lambda2);
                targetFunctionCalled += 6;
            }
            else
            {
                rightBorder = lambda2;
                lambda2 = lambda1;
                f2 = f1;
                interval = rightBorder - leftBorder;
                lambda1 = rightBorder - T * interval;
                f1 = func(lambda1);
                targetFunctionCalled += 6;
            }
            cycles++;
            targetFunctionCalled++;
        }

        return (f1 < f2) ? lambda1 : lambda2;
    }

    public static (double, double, double) PartialGradientsWithPenalty(double x, double y, double z, double r)
    {
        double gradX_f = y * z;
        double gradY_f = x * z;
        double gradZ_f = x * y;

        //Console.WriteLine($"{gradX_f}, {gradY_f}, {gradZ_f}");

        double gradX_g = 2 * (y + z);
        double gradY_g = 2 * (x + z);
        double gradZ_g = 2 * (x + y);
        //Console.WriteLine($"{gradX_g}, {gradY_g}, {gradZ_g}");

        double gradX = gradX_f + 2 * r * EqualityConstraintValue(x, y, z) * gradX_g;
        double gradY = gradY_f + 2 * r * EqualityConstraintValue(x, y, z) * gradY_g;
        double gradZ = gradZ_f + 2 * r * EqualityConstraintValue(x, y, z) * gradZ_g;
        //Console.WriteLine($"{gradZ}, {gradY}, {gradZ}");

        return (gradX, gradY, gradZ);
    }
    public static (double, double, double) PartialGradients(double x, double y, double z)
    {
        return (y * z, x * z, x * y);
    }

    public static double PenaltyFunction(double x, double y, double z, double r)
    {
        double target = TargetFunction(x, y, z);
        double equalityPenalty = Math.Pow(EqualityConstraintValue(x, y, z), 2);

        // might be incorrect due to MIN(0, x) instead of MAX(0, x) or even it should be a different method
        double inequalityPenalty = Math.Pow(Math.Min(0, x), 2) + Math.Pow(Math.Min(0, y), 2) + Math.Pow(Math.Min(0, z), 2); 
        //double inequalityPenalty = Pow(Min(x, 0),2) + Pow(Min(y, 0), 2) + Pow(Min(z, 0), 2);

        return target + r * (equalityPenalty + inequalityPenalty);
    }

    public static double TargetFunction(double x, double y, double z)
    {
        return -x * y * z;
    }

    public static double EqualityConstraintValue(double x, double y, double z)
    {
        //x = Math.Abs(x);
        //y = Math.Abs(y);
        //z = Math.Abs(z);
        return (x * y + y * z + x * z) * 2 - 1;
        //return Math.Abs((x*y + y*z + x*z) *2 -1) < Epsilon;
    }

    public static (bool, double, double, double) NegativeConstarint(double x, double y, double z)
    {
        double nx = x > 0 ? x : 0;
        double ny = y > 0 ? y : 0;
        double nz = z > 0 ? z : 0;
        bool isConstraintViolated = x >= 0 && y >= 0 && z >= 0;

        return (isConstraintViolated, nx, ny, nz);
    }
}

