using Microsoft.VisualBasic;
using System.Drawing;
using System.Numerics;
using Plotly.NET;
using Plotly.NET.LayoutObjects;
using System;
using System.Collections.Generic;
using System.Linq;
using static Plotly.NET.StyleParam;


class Program()
{
    public const double Epsilon = 1e-2;
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
        //z = 0;

        //double r = 1;
        //Func<double, double, double, double, double> targetFunctionFunc = (x, y, z, r) => TargetFunctionWithPenalty(x, y, z, r);

        ////var (pointTimeline, functionCalls, cycles) = FastestDescentWithPenalty(x, y, z, r, targetFunctionFunc);
        //var (pointTimeline, functionCalls, cycles) = GradientDescentWithPenalty(x, y, z, r, targetFunctionFunc);

        //Console.WriteLine("==========================================");
        //int count = 0;
        //foreach (var point in pointTimeline)
        //{
        //    Console.WriteLine($"\n{count} ({point.Item1}, {point.Item2}, {point.Item3}), lr: {point.Item4}, f: {point.Item5}");
        //    count++;
        //}
        //Console.WriteLine($"Function Calls: {functionCalls}, Cycles: {cycles}");

        TranslateFromPython.Mainpiton(x,y,z,1,5);
    }

    public static void GradientDescentWithPenalty(
        double x, double y, double z, double r, Func<double, double, double, double, double> targetFunctionFunc)
    {
        List<(double, double, double, double, double)> pointTimeLine = new();
        pointTimeLine.Add((x, y, z, 1, targetFunctionFunc(x, y, z, r)));
        int cycles = 0;
        int iterations = 0;
        int targetFunctionCalled = 0;
        double lr = 0.01;
        double old_b = PenaltyFunction(x,y,z,r);

        while (true)
        {

            (double gradX, double gradY, double gradZ) = PartialGradientsWithPenalty(x,y,z,r);

            double b = PenaltyFunction(x,y,z,r);
            targetFunctionCalled++;

            x -= lr * gradX;
            y -= lr * gradY;
            z -= lr * gradZ;
            Console.WriteLine($"{iterations}| ({x}, {y}, {z}) | {b}");
            r /= 2;
            lr /= 2;
            cycles++;

            if (old_b < b)
            {
                iterations++;
                if(iterations >= 5)
                {
                    Console.WriteLine($"Mehtod stoped due to b getting worse");
                    break;
                }
            }
            old_b = b;
        }

        //var last = pointTimeLine.Last();
        //Console.WriteLine($"{last.Item1}, {last.Item2}, {last.Item3}, {targetFunctionCalled * 6}, {cycles}");
    }

    public static (List<(double, double, double, double, double)>, int, int) FastestDescentWithPenalty(
        double x, double y, double z, double r, Func<double, double, double, double, double> targetFunctionFunc)
    {
        List<(double, double, double, double, double)> pointTimeLine = new();
        pointTimeLine.Add((x, y, z, 1, targetFunctionFunc(x, y, z, r)));
        int cycles = 0;
        int iterations = 0;
        int targetFunctionCalled = 0;
        double old_f = PenaltyFunction(x, y, z, r);

        while (iterations < 20)
        {
            (double gradX, double gradY, double gradZ) = PartialGradientsWithPenalty(x, y, z, r);
            targetFunctionCalled++;
            Func<double, double> targetFunctionValueFunc = lambda => targetFunctionFunc(x - lambda * gradX, y - lambda * gradY, z - lambda * gradZ, r);

            double lambda = GoldenSectionSearch(targetFunctionValueFunc, ref cycles, ref targetFunctionCalled);
            //double lambda = HalfSectionSearch(targetFunctionValueFunc, ref cycles, ref targetFunctionCalled);

            //Console.WriteLine(x < lambda * gradX);
            x -= lambda * gradX;
            y -= lambda * gradY;
            z -= lambda * gradZ;
            r /= 5;
            cycles++;
            iterations++;

            double new_f = PenaltyFunction(x, y, z, r);

            if (iterations == 30){
                break;
            }
            //if (old_f < new_f)
            //{
            //    Console.WriteLine($"\n\nold_f < new_f: {old_f < new_f} | {old_f} < {new_f}");
            //    Console.WriteLine($"{iterations}:, {x}, {y}, {z}, r: {r} | {new_f} , {TargetFunctionWithPenalty(x,y,z,r)}\n");
            //    break;
            //}

            //if (Math.Abs(EqualityConstraintValue(x,y,z)) < Epsilon) {
            //    Console.WriteLine($"\n Eq= {EqualityConstraintValue(x,y,z)} maziau nei Epsilon\n");
            //    break;
            //}
            old_f = new_f;
            Console.WriteLine($"\n{iterations}:, {x}, {y}, {z}, r: {r} | {new_f}");
            pointTimeLine.Add((x, y, z, lambda, TargetFunctionNotMinimized(x, y, z)));

        }

        return (pointTimeLine, targetFunctionCalled * 6, cycles);
    }


    private static double HalfSectionSearch(Func<double, double> targetFunction,
        ref int cycles, ref int targetFunctionCalled, double leftBorder = 0, double rightBorder = 2)
    {
        double interval = rightBorder - leftBorder;
        double middlePoint = (leftBorder + rightBorder) / 2;
        double middlePointValue = targetFunction(middlePoint);
        double lambda = middlePoint;

        while (interval > Epsilon)
        {
            double leftPoint = leftBorder + interval / 4;
            double leftPointValue = targetFunction(leftPoint);
            targetFunctionCalled++;
            if (leftPointValue < middlePointValue)
            {
                rightBorder = middlePoint;
                middlePoint = leftPoint;
                middlePointValue = leftPointValue;
            }
            else
            {
                double rightPoint = rightBorder - interval / 4;
                double rightPointValue = targetFunction(rightPoint);
                targetFunctionCalled++;
                leftBorder = middlePoint;
                middlePoint = rightPoint;
                middlePointValue = rightPointValue;
            }
            interval = rightBorder - leftBorder;
            lambda = middlePoint;
            cycles++;
        }
        return lambda;
    }

    private static double GoldenSectionSearch(Func<double, double> func,
        ref int cycles, ref int targetFunctionCalled, double leftBorder = 0, double rightBorder = 4)
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

    //public static (double, double, double) PartialGradientsWithPenalty(double x, double y, double z, double r)
    //{
    //    (double gradientX, double gradientY, double gradientZ) = PartialGradients(x, y, z);

    //    double equalityConstraint = EqualityConstraintValue(x, y, z);
    //    double ecX = equalityConstraint * 4 * (y + z);
    //    double ecY = equalityConstraint * 4 * (x + z);
    //    double ecZ = equalityConstraint * 4 * (x + y);

    //    double icX = x < 0 ? (2 * x) : 0;
    //    double icY = y < 0 ? (2 * y) : 0;
    //    double icZ = z < 0 ? (2 * z) : 0;

    //    gradientX += (icX + ecX) / r;
    //    gradientY += (icY + ecY) / r;
    //    gradientZ += (icZ + ecZ) / r;

    //    return (gradientX, gradientY, gradientZ);
    //}
    public static (double, double, double) PartialGradientsWithPenalty(double x, double y, double z, double r)
    {
        (double gradientX, double gradientY, double gradientZ) = PartialGradients(x, y, z);

        double equalityConstraint = EqualityConstraintValue(x, y, z);
        double ecX = 2 * Math.Abs(y + z);
        double ecY = 2 * Math.Abs(x + z);
        double ecZ = 2 * Math.Abs(x + y);

        double icX = x < 0 ? x : 0;
        double icY = y < 0 ? y : 0;
        double icZ = z < 0 ? z : 0;

        gradientX += (icX + ecX) / r;
        gradientY += (icY + ecY) / r;
        gradientZ += (icZ + ecZ) / r;

        return (gradientX, gradientY, gradientZ);
    }

    public static double TargetFunctionWithPenalty(double x, double y, double z, double r)
    {
        return TargetFunction(x, y, z) + PenaltyFunction(x, y, z, r);
    }

    //private static double PenaltyFunction(double x, double y, double z, double r)
    //{
    //    double equalityConstraint = Math.Pow(EqualityConstraintValue(x, y, z), 2);
    //    Func<double, double> inequalityFunc = (x) => (Math.Pow(AvoidNegative(x), 2));
    //    double inequalityConstraint = inequalityFunc(x) + inequalityFunc(y) + inequalityFunc(z);

    //    return (1 / r) * (equalityConstraint + inequalityConstraint);
    //}
    private static double PenaltyFunction(double x, double y, double z, double r)
    {
        double equalityConstraint = Math.Abs(EqualityConstraintValue(x, y, z));
        Func<double, double> inequalityFunc = (x) => AvoidNegative(x);
        double inequalityConstraint = inequalityFunc(x) + inequalityFunc(y) + inequalityFunc(z);

        return (1 / r) * (equalityConstraint + inequalityConstraint);
    }

    public static double EqualityConstraintValue(double x, double y, double z)
    {
        return (x * y + y * z + x * z) * 2 - 1;
    }

    public static double AvoidNegative(double x)
    {
        //Console.WriteLine($"{x}| {Math.Abs(Math.Min(x, 0))}, {Math.Min(x, 0)}, {Math.Max(-x, 0)}");
        //return Math.Abs(Math.Min(x, 0));
        //return Math.Min(x, 0);
        return Math.Max(-x, 0);
    }

    public static double TargetFunction(double x, double y, double z)
    {
        return -x * y * z;
    }

    public static double TargetFunctionNotMinimized(double x, double y, double z)
    {
        return x * y * z;
    }

    public static (double, double, double) PartialGradients(double x, double y, double z)
    {
        return (-y * z, -x * z, -x * y);
    }

    //public static (List<(double, double, double, double, double)>, int, int) GradientDescentWithPenalty(
    //    double x, double y, double z, double r, Func<double, double, double, double, double> targetFunctionFunc)
    //{
    //    List<(double, double, double, double, double)> pointTimeLine = new();
    //    pointTimeLine.Add((x, y, z, 1, targetFunctionFunc(x, y, z, r)));
    //    int cycles = 0;
    //    int iterations = 0;
    //    int targetFunctionCalled = 0;
    //    double old_f = PenaltyFunction(x, y, z, r);
    //    double learningRate = 0.1;

    //    while (true)
    //    {
    //        double new_f = PenaltyFunction(x, y, z, r);
    //        if (new_f > old_f)
    //        {
    //            Console.WriteLine($"old_f < new_f: {old_f} < {new_f}");
    //        }
    //        old_f = new_f;
    //        //if (!(r > Epsilon || x < 0 || y < 0 || z < 0))
    //        //{
    //        //    Console.WriteLine($"Stopped. r: {r}, x: {x}, y: {y}, z: {z}");
    //        //    break;
    //        //}
    //        (double gradX, double gradY, double gradZ) = PartialGradientsWithPenalty(x, y, z, r);
    //        targetFunctionCalled++;

    //        x -= learningRate * gradX;
    //        y -= learningRate * gradY;
    //        z -= learningRate * gradZ;

    //        pointTimeLine.Add((x, y, z, 1, targetFunctionFunc(x, y, z, r)));
    //        Console.WriteLine($"Iteration: {iterations}, x: {x}, y: {y}, z: {z}, r: {r}\n");

    //        if (iterations > 5)
    //        {
    //            r = Math.Max(r / 5, Epsilon);
    //            learningRate /= 5;
    //        }

    //        if (Math.Sqrt(gradX * gradX + gradY * gradY + gradZ * gradZ) < Epsilon)
    //        {
    //            break;
    //        }

    //        cycles++;
    //        iterations++;
    //    }

    //    return (pointTimeLine, targetFunctionCalled * 3, cycles);
    //}
}

