using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Security.Cryptography.X509Certificates;

class Program() {
    public const double Epsilon = 1e-6;
    public static void Main(string[] args)
    {
        double x = 0;
        double y = 0;
        double z = 0;

        //double x = 1;
        //double y = 1;
        //double z = 1;

        //double x = 0.5;
        //double y = 0.8;
        //double z = 0.0;



    }

    public static (List<(double, double, double, double)>, int, int) FastestDescent(double x, double y, double z)
    {
        List<(double, double, double, double)> pointTimeLine = new();
        pointTimeLine.Add((x, y, z, 1));
        int cycles = 0;
        int targetFunctionCalled = 0;

        while (cycles < 50000)
        {
            (double gradX, double gradY, double gradZ) = PartialGradients(x, y, z);
            targetFunctionCalled += 3;
            if (Math.Abs(gradX) < Epsilon && Math.Abs(gradY) < Epsilon && Math.Abs(gradZ) < Epsilon)
            {
                Console.WriteLine($"Metodas pasibaigė, nes gradientų reikšmės mažesnės" +
                    $" už Epsilon {Epsilon}: gradX: {gradX}, gradY: {gradY}, gradZ: {gradZ}\n");
                break;
            }

            Func<double, double> goldenSectionFunc = lambda => TargetFunction(x - lambda * gradX, y - lambda * gradY, z - lambda * gradZ);

            double lambda= GoldenSectionSearch(goldenSectionFunc, ref cycles, ref targetFunctionCalled);
            x -= lambda * gradX;
            y -= lambda * gradY;
            z -= lambda * gradZ;

            pointTimeLine.Add((x, y, z, lambda));
            cycles++;
        }

        return (pointTimeLine, targetFunctionCalled, cycles);
    }

    private static double Penalty(double x, double y, double z, double penalty)
    {
        return TargetFunction(x, y, z) + penalty * Math.Pow(EqualityConstraint(x, y, z) ? 0 : (x * y + y * z + x * z) * 2 - 1, 2);
    }

    private static double GoldenSectionSearch(Func<double, double> func, ref int cycles, ref int targetFunctionCalled, double leftBorder = 0, double rightBorder = 5)
    {
        double T = (-1 + Math.Sqrt(5)) / 2;
        double interval = rightBorder - leftBorder;
        double lambda1 = rightBorder - T * interval;
        double lambda2 = leftBorder + T * interval;

        double f1 = func(lambda1);
        double f2 = func(lambda2);
        targetFunctionCalled += 2;

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
            }
            else
            {
                rightBorder = lambda2;
                lambda2 = lambda1;
                f2 = f1;
                interval = rightBorder - leftBorder;
                lambda1 = rightBorder - T * interval;
                f1 = func(lambda1);
            }
            cycles++;
            targetFunctionCalled++;
        }

        return (f1 < f2) ? lambda1 : lambda2;
    }


    public static double TargetFunction(double x, double y, double z)
    {
        return -x*y*z;
    }

    public static (double, double, double) PartialGradients(double x, double y, double z)
    {
        return (y*z, x*z, x*y);
    }

    public static bool EqualityConstraint(double x, double y, double z)
    {
        return Math.Abs((x*y + y*z + x*z) *2 -1) < Epsilon;
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

