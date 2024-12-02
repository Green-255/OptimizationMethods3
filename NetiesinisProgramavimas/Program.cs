class Program()
{
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

        int maxIterations = 1000;
        double reductionFactor = 0.1;
        double penaltyMultiplier = 10.0;

        for (int i = 0; i < maxIterations; i++)
        {
            Console.WriteLine($"Iteration {i + 1}, Penalty Multiplier: {penaltyMultiplier}");
            Func<double, double, double, double, double> penaltyFunc = (px, py, pz, pr) => PenaltyFunction(px, py, pz, pr);

            var (pointTimeline, functionCalls, cycles) = FastestDescentWithPenalty(x, y, z, penaltyMultiplier, penaltyFunc);
            var (finalX, finalY, finalZ, _) = pointTimeline[^1];

            Console.WriteLine($"Final Point: x = {finalX}, y = {finalY}, z = {finalZ}");
            Console.WriteLine($"Function Calls: {functionCalls}, Cycles: {cycles}");

            x = finalX;
            y = finalY;
            z = finalZ;

            penaltyMultiplier *= reductionFactor;
        }
    }

    public static (List<(double, double, double, double)>, int, int) FastestDescentWithPenalty(
        double x, double y, double z, double penaltyMultiplier, Func<double, double, double, double, double> penaltyFunc)
    {
        List<(double, double, double, double)> pointTimeLine = new();
        pointTimeLine.Add((x, y, z, 1));
        int cycles = 0;
        int targetFunctionCalled = 0;

        while (cycles < 50000)
        {
            (double gradX, double gradY, double gradZ) = PartialGradientsWithPenalty(x, y, z, penaltyMultiplier, penaltyFunc);
            targetFunctionCalled += 3;

            if (Math.Abs(gradX) < Epsilon && Math.Abs(gradY) < Epsilon && Math.Abs(gradZ) < Epsilon)
            {
                Console.WriteLine($"Converged: gradX: {gradX}, gradY: {gradY}, gradZ: {gradZ}");
                break;
            }

            Func<double, double> goldenSectionFunc = lambda => penaltyFunc(x - lambda * gradX, y - lambda * gradY, z - lambda * gradZ, penaltyMultiplier);

            double lambda = GoldenSectionSearch(goldenSectionFunc, ref cycles, ref targetFunctionCalled);
            x -= lambda * gradX;
            y -= lambda * gradY;
            z -= lambda * gradZ;

            pointTimeLine.Add((x, y, z, lambda));
            cycles++;
        }

        return (pointTimeLine, targetFunctionCalled, cycles);
    }

    public static double PenaltyFunction(double x, double y, double z, double r)
    {
        double target = TargetFunction(x, y, z);
        double equalityPenalty = Math.Pow(EqualityConstraint(x, y, z), 2);
        double inequalityPenalty = Math.Pow(Math.Min(0, x), 2) + Math.Pow(Math.Min(0, y), 2) + Math.Pow(Math.Min(0, z), 2);

        return target + r * (equalityPenalty + inequalityPenalty);
    }


    public static (List<(double, double, double, double)>, int, int) FastestDescent(
        double x, double y, double z, double penaltyMultiplier, Func<double, double, double, double, double> penaltyFunc)
    {
        List<(double, double, double, double)> pointTimeLine = new();
        pointTimeLine.Add((x, y, z, 1));
        int cycles = 0;
        int targetFunctionCalled = 0;

        while (cycles < 50000)
        {
            (double gradX, double gradY, double gradZ) = PartialGradientsWithPenalty(x, y, z, penaltyMultiplier, penaltyFunc);
            targetFunctionCalled += 3;

            if (Math.Abs(gradX) < Epsilon && Math.Abs(gradY) < Epsilon && Math.Abs(gradZ) < Epsilon)
            {
                Console.WriteLine($"Converged: gradX: {gradX}, gradY: {gradY}, gradZ: {gradZ}");
                break;
            }

            Func<double, double> goldenSectionFunc = lambda => penaltyFunc(x - lambda * gradX, y - lambda * gradY, z - lambda * gradZ, penaltyMultiplier);

            double lambda = GoldenSectionSearch(goldenSectionFunc, ref cycles, ref targetFunctionCalled);
            x -= lambda * gradX;
            y -= lambda * gradY;
            z -= lambda * gradZ;

            pointTimeLine.Add((x, y, z, lambda));
            cycles++;
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

    public static (double, double, double) PartialGradientsWithPenalty(double x, double y, double z, double r, Func<double, double, double, double, double> penaltyFunc)
    {
        double delta = 1e-5;
        double gradX = (penaltyFunc(x + delta, y, z, r) - penaltyFunc(x - delta, y, z, r)) / (2 * delta);
        double gradY = (penaltyFunc(x, y + delta, z, r) - penaltyFunc(x, y - delta, z, r)) / (2 * delta);
        double gradZ = (penaltyFunc(x, y, z + delta, r) - penaltyFunc(x, y, z - delta, r)) / (2 * delta);

        return (gradX, gradY, gradZ);
    }
    public static (double, double, double) PartialGradients(double x, double y, double z)
    {
        return (y * z, x * z, x * y);
    }


    public static double TargetFunction(double x, double y, double z)
    {
        return -x * y * z;
    }

    public static double EqualityConstraint(double x, double y, double z)
    {
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

