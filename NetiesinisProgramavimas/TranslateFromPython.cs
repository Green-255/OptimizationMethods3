using System;
using System.Collections.Generic;
using System.Linq;

class TranslateFromPython
{
    public static void Mainpiton(double a, double b, double c, double r, double rChange)
    {
        var (x, y, z, cycles, functionCalls) = OptimizeFunction((a, b, c), r, rChange);
        Console.WriteLine($"Result: x = {x}, y = {y}, z = {z}\nFunkcijos kvietimų: {functionCalls}, Ciklų kiekis: {cycles}");
        Console.WriteLine($"{TargetFunctionNotMinimized(x, y, z)}");
    }

    public static double Distance((double, double, double) point1, (double, double, double) point2)
    {
        return Math.Sqrt(
            Math.Pow(point2.Item1 - point1.Item1, 2) +
            Math.Pow(point2.Item2 - point1.Item2, 2) +
            Math.Pow(point2.Item3 - point1.Item3, 2)
        );
    }

    public static (double, double, double, int, int) OptimizeFunction((double, double, double) start, double r, double rChange)
    {
        double x = start.Item1;
        double y = start.Item2;
        double z = start.Item3;

        double lr = 0.01;
        double epsilon = 1e-6;
        int functionCalls = 0;
        int cycles = 0;
        int i = 0;
        double old_b = Int32.MaxValue;

        while (true)
        {
            Console.WriteLine($"{x} & {y} & {z} \\\\");
            //Console.WriteLine($"{r} & {old_b} \\\\");

            (x, y, z, int cyclesAddition, int functionCallsAddition) = GradientDescent((x, y, z), lr, epsilon, r);
            functionCalls += functionCallsAddition;
            cycles += cyclesAddition;
            i++;

            r /= rChange;
            epsilon /= rChange;
            lr /= rChange;

            double new_b = PenaltyFunction(x, y, z, r);
            functionCalls += 2;
            if (old_b <= new_b)
            {
                break;
            }
            old_b = new_b;
        }

        return (x, y, z, cycles, functionCalls);
    }

    public static (double, double, double, int, int) GradientDescent((double, double, double) X, double lr, double epsilon, double r)
    {
        int cycles = 0;
        int functionCalls = 0;

        var (x, y, z) = (X.Item1, X.Item2, X.Item3);
        var (newX, newY, newZ) = (x, y, z);

        do
        {
            x = newX;
            y = newY;
            z = newZ;

            (double gradX, double gradY, double gradZ) = PartialGradientsWithPenalty(x, y, z, r);
            newX = x - lr * gradX;
            newY = y - lr * gradY;
            newZ = z - lr * gradZ;

            functionCalls += 8;
            cycles++;
        } while (Distance((x, y, z), (newX, newY, newZ)) >= epsilon);

        return (newX, newY, newZ, cycles, functionCalls);
    }

    private static double PenaltyFunction(double x, double y, double z, double r)
    {
        Func<double, double> inequalityFunc = (x) => (Math.Pow(NegativePenaklty(x), 2));
        double equalityConstraint = Math.Pow(EqualityConstraintValue(x, y, z), 2);
        double inequalityConstraint = inequalityFunc(x) + inequalityFunc(y) + inequalityFunc(z);

        return (1 / r) * (inequalityConstraint + equalityConstraint);
    }


    public static (double, double, double) PartialGradientsWithPenalty(double x, double y, double z, double r)
    {
        double equalityConstraint = EqualityConstraintValue(x, y, z);
        double ecX = equalityConstraint * 4 * (y + z);
        double ecY = equalityConstraint * 4 * (x + z);
        double ecZ = equalityConstraint * 4 * (x + y);

        double icX = x < 0 ? (2 * x) : 0;
        double icY = y < 0 ? (2 * y) : 0;
        double icZ = z < 0 ? (2 * z) : 0;

        (double gradientX, double gradientY, double gradientZ) = PartialGradients(x, y, z);

        gradientX += (icX + ecX) / r;
        gradientY += (icY + ecY) / r;
        gradientZ += (icZ + ecZ) / r;

        return (gradientX, gradientY, gradientZ);
    }

    public static double EqualityConstraintValue(double x, double y, double z)
    {
        return (x * y + y * z + x * z) * 2 - 1;
    }

    public static (double, double, double) PartialGradients(double x, double y, double z)
    {
        return (-y * z, -x * z, -x * y);
    }
    public static double NegativePenaklty(double x)
    {
        return Math.Max(-x, 0);
    }
    public static double TargetFunctionNotMinimized(double x, double y, double z)
    {
        return x * y * z;
    }
}
