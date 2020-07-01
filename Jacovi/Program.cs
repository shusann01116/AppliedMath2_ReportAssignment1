using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Jacovi {
    class Program {

        static void Main(string[] args) {
            var loadDataCsv = @"../../../Resource/data.csv";
            var saveDataCsv = @"../../../Resource/savedata.csv";

            var data = LoadCSV(loadDataCsv);

            var dataset = new List<List<int[]>> {
                new List<int[]>(),
                new List<int[]>(),
            };

            data = data.OrderBy(x => Guid.NewGuid()).ToList();
            dataset[0] = data.GetRange(0, data.Count / 2);
            dataset[1] = data.GetRange((int)Math.Ceiling(data.Count / 2.0), data.Count / 2);

            var returnValues = new ReturnValues();
            foreach (var m in Enumerable.Range(3, 5)) {
                Console.WriteLine(m);
                returnValues.Values.Add(Executer(m, data, data));
                returnValues.Values.Add(Executer(m, dataset[0], dataset[1]));
                returnValues.Values.Add(Executer(m, dataset[1], dataset[0]));
            }

            returnValues.SaveToCSV(saveDataCsv);

            Console.ReadKey();
        }

        static ReturnValue Executer(int m, List<int[]> data, List<int[]> validator) {
            var A = SolveMatrixEquation(m, data);
            var error = CalculateError(data, A);
            var validation = CalculateError(validator, A);
            var penalty = CalculatePenalty(A);

            Console.WriteLine($"Error:{error}, Validation:{validation}, Cross ratio:{error / validation}, Penalty:{penalty}");
            return new ReturnValue {M = m, A = A, Error = error, Penalty = penalty, Validation = validation};
        }



        static double[,] SolveMatrixEquation(int m, List<int[]> data) {
            // データ数
            var n = data.Count;

            // X行列
            var X = new double[m, m];

            // Yベクトル
            var Y = new double[m, 1];

            // X_jkを導出
            for (var j = 1; j <= m; j++) {
                for (var k = 1; k <= m; k++) {
                    double summation = 0;
                    for (var i = 1; i <= n; i++) {
                        var x = data[i - 1][0];
                        summation += Math.Pow(x, 2 * m - j - k);
                    }

                    X[j - 1, k - 1] = summation;
                }
            }

            // Y_jを導出
            for (var j = 1; j <= m; j++) {
                double summation = 0;
                for (var i = 1; i <= n; i++) {
                    var x = data[i - 1][0];
                    var y = data[i - 1][1];
                    summation += Math.Pow(x, m - j) * y;
                }

                Y[j - 1, 0] = summation;

            }

            var XInverse = InverseMatrix(X);
            return MatrixTimesMatrix(XInverse, Y);
        }

        static double CalculateError(List<int[]> data, double[,] A) {
            double error = 0;

            double CalculateModelEquation(int x) {
                double result = 0;
                for (int i = 1; i <= A.Length; i++) {
                    result += A[i - 1, 0] * Math.Pow(x, A.Length - i);
                }

                return result;
            }

            for (int i = 0; i < A.Length; i++) {
                error += Math.Pow(data[i][1] - CalculateModelEquation(data[i][0]), 2);
            }

            return error;
        }

        static double CalculatePenalty(double[,] A) {
            double result = 0;
            for (int i = 0; i < A.Length; i++) {
                result += Math.Pow(A[i, 0], 2);
            }

            return result;
        }

        static List<int[]> LoadCSV(string filePath) {
            var data = new List<int[]>();
            var fs = File.OpenRead(filePath);
            using var sr = new StreamReader(fs);
            while (!sr.EndOfStream) {
                var line = sr.ReadLine().Split(',');
                data.Add(new[] {
                    int.Parse(line[0]),
                    int.Parse(line[1]),
                });
            }

            return data;
        }

        static double[,] InverseMatrix(double[,] A) {

            int n = A.GetLength(0);
            int m = A.GetLength(1);

            double[,] invA = new double[n, m];

            if (n == m) {

                int max;
                double tmp;

                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < n; i++) {
                        invA[j, i] = (i == j) ? 1 : 0;
                    }
                }

                for (int k = 0; k < n; k++) {
                    max = k;
                    for (int j = k + 1; j < n; j++) {
                        if (Math.Abs(A[j, k]) > Math.Abs(A[max, k])) {
                            max = j;
                        }
                    }

                    if (max != k) {
                        for (int i = 0; i < n; i++) {
                            // 入力行列側
                            tmp = A[max, i];
                            A[max, i] = A[k, i];
                            A[k, i] = tmp;
                            // 単位行列側
                            tmp = invA[max, i];
                            invA[max, i] = invA[k, i];
                            invA[k, i] = tmp;
                        }
                    }

                    tmp = A[k, k];

                    for (int i = 0; i < n; i++) {
                        A[k, i] /= tmp;
                        invA[k, i] /= tmp;
                    }

                    for (int j = 0; j < n; j++) {
                        if (j != k) {
                            tmp = A[j, k] / A[k, k];
                            for (int i = 0; i < n; i++) {
                                A[j, i] = A[j, i] - A[k, i] * tmp;
                                invA[j, i] = invA[j, i] - invA[k, i] * tmp;
                            }
                        }
                    }

                }


                //逆行列が計算できなかった時の措置
                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < n; i++) {
                        if (double.IsNaN(invA[j, i])) {
                            Console.WriteLine("Error : Unable to compute inverse matrix");
                            invA[j, i] = 0;//ここでは，とりあえずゼロに置き換えることにする
                        }
                    }
                }


                return invA;

            } else {
                Console.WriteLine("Error : It is not a square matrix");
                return invA;
            }

        }

        static double[,] MatrixTimesMatrix(double[,] A, double[,] B) {

            double[,] product = new double[A.GetLength(0), B.GetLength(1)];

            for (int i = 0; i < A.GetLength(0); i++) {
                for (int j = 0; j < B.GetLength(1); j++) {
                    for (int k = 0; k < A.GetLength(1); k++) {
                        product[i, j] += A[i, k] * B[k, j];
                    }
                }
            }

            return product;

        }

        #region Jacovi

        static void Jacovi() {
            var filePath = @"../../../Resource/data.csv";

            var data = LoadCSV(filePath);

            var m = 3;
            var A = new double[m, 1];

            var result = InitializeJacovi(m, data);

            foreach (var i in Enumerable.Range(0, 100)) {
                A = ExecuteJacovi(m, result, A);
                if (i % 5 == 0) {
                    Console.WriteLine($"{A[0, 0]:G}, {A[1, 0]:G}, {A[2, 0]:G}");
                }
            }
        }

        /// <summary>
        /// ヤコビ法で未定係数を決定する
        /// </summary>
        /// <param name="m">モデル関数の項数</param>
        static List<double[,]> InitializeJacovi(int m, List<int[]> data) {
            // データ数
            var n = data.Count;

            // X行列
            var X = new double[m, m];

            // Yベクトル
            var Y = new double[m, 1];

            // X_jkを導出
            for (var j = 1; j <= m; j++) {
                for (var k = 1; k <= m; k++) {
                    double summation = 0;
                    for (var i = 1; i <= n; i++) {
                        var x = data[i - 1][0];
                        summation += Math.Pow(x, 2 * m - j - k);
                    }

                    X[j - 1, k - 1] = summation;
                }
            }

            // Y_jを導出
            for (var j = 1; j <= m; j++) {
                double summation = 0;
                for (var i = 1; i <= n; i++) {
                    var x = data[i - 1][0];
                    var y = data[i - 1][1];
                    summation += Math.Pow(x, m - j) * y;
                }

                Y[j - 1, 0] = summation;
            }

            // 対角行列D
            var D = new double[m, m];
            for (var i = 0; i < m; i++) {
                D[i, i] = X[i, i];
            }

            // 上三角行列U+下三角行列L
            var UplusL = new double[m, m];

            #region UplusL = X

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    UplusL[i, j] = X[i, j];
                }
            }

            #endregion
            for (var i = 0; i < m; i++) {
                UplusL[i, i] = 0;
            }

            var D_i = InverseMatrix(D);
            // D^-1*Y
            var D_iTimesY = MatrixTimesMatrix(D_i, Y);
            // D^-1*(U+L)
            var D_iTimesUplusL = MatrixTimesMatrix(D_i, UplusL);

            var result = new List<double[,]> {
                D_iTimesY,
                D_iTimesUplusL
            };

            return result;
        }

        /// <summary>
        /// InitilizeJacoviの結果を用いて未定係数を決定する
        /// </summary>
        /// <param name="m">モデル関数の項数</param>
        /// <param name="InitializedResult">InitilizeJacoviの結果を代入</param>
        /// <param name="A">未定係数(m行1列のベクトル)</param>
        /// <returns></returns>
        static double[,] ExecuteJacovi(int m, List<double[,]> InitializedResult, double[,] A) {
            var A_next = new double[m, 1];
            var D_iTimesY = InitializedResult[0];
            var D_iTimesUplusL = InitializedResult[1];

            // D^-1*(U+L)*a^q
            var D_iTimesUplusLTimesA = MatrixTimesMatrix(D_iTimesUplusL, A);
            for (var i = 0; i < m; i++) {
                A_next[i, 0] = D_iTimesY[i, 0] - D_iTimesUplusLTimesA[i, 0];
            }

            return A_next;
        }

        #endregion
    }

    class ReturnValue {
        public int M { get; set; }
        public double[,] A { get; set; }
        public double Error { get; set; }
        public double Validation { get; set; }
        public double Penalty { get; set; }
    }

    class ReturnValues {
        public IList<ReturnValue> Values { get; set; } = new List<ReturnValue>();
        public void SaveToCSV(string filePath) {
            using var fs = File.OpenWrite(filePath);
            using var sr = new StreamWriter(fs);
            sr.WriteLine(@"M, Error, Validation, Cross ratio, Penalty, A");

            var sb = new StringBuilder();
            foreach (var v in Values) {
                sb.Append($"{v.M}, {v.Error}, {v.Validation}, {v.Validation / v.Error}, {v.Penalty}");
                foreach (var d in v.A) {
                    sb.Append($", {d}");
                }

                sr.WriteLine(sb);
                sb.Clear();
            }
        }
    }
}
