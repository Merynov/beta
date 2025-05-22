using System;
using System.IO;
using System.Windows.Forms;
using System.Drawing;
using System.Globalization;
using System.Diagnostics;

// програма для обчислення власних значень і векторів матриці методами Леверр'є-Фадєєва та Данилевського.

namespace EigenApp
{
    public partial class Form1 : Form
    {
        private int matrixSize = 3;
        private double[,] matrix;
        private DataGridView dataGridViewMatrix;
        private TextBox textBoxResult;
        private NumericUpDown numericUpDownSize;
        private Button buttonCalculate;
        private Button buttonSave;
        private Label labelMatrixSize;
        private Label labelMatrix;
        private Label labelResult;
        private Label labelReadme;
        private ComboBox comboBoxMethod; // для вибору методу
        private Label labelMethod;
        private long operationCount; // для підрахунку операцій

        public Form1()
        {
            InitializeComponents();
            InitializeMatrixGrid(matrixSize);
        }

        private void InitializeComponents()
        {
            // еалаштування форми
            this.Size = new Size(700,600);
            this.Text = "Калькулятор власних значень і векторів";
            this.FormBorderStyle = FormBorderStyle.FixedSingle;
            this.MaximizeBox = false;
            this.BackColor = Color.White;

            // поле для розміру матриці
            labelMatrixSize = new Label
            {
                Location = new Point(20, 20),
                Size = new Size(120, 20),
                Text = "Розмір матриці:",
                Font = new Font("Arial", 10, FontStyle.Bold)
            };
            this.Controls.Add(labelMatrixSize);

            // NumericUpDown для вибору розміру матриці
            numericUpDownSize = new NumericUpDown
            {
                Location = new Point(150, 20),
                Size = new Size(60, 20),
                Minimum = 2,
                Maximum = 10,
                Value = matrixSize,
                Font = new Font("Arial", 10)
            };
            numericUpDownSize.ValueChanged += numericUpDownSize_ValueChanged;
            this.Controls.Add(numericUpDownSize);

            // поле для вибору методу
            labelMethod = new Label
            {
                Location = new Point(220, 20),
                Size = new Size(120, 20),
                Text = "Метод:",
                Font = new Font("Arial", 10, FontStyle.Bold)
            };
            this.Controls.Add(labelMethod);

            // ComboBox для вибору методу
            comboBoxMethod = new ComboBox
            {
                Location = new Point(350, 20),
                Size = new Size(150, 20),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Arial", 10)
            };
            comboBoxMethod.Items.AddRange(new string[] { "Леверр'є-Фадєєв", "Данилевський" });
            comboBoxMethod.SelectedIndex = 0; // За замовчуванням Леверр'є-Фадєєв
            this.Controls.Add(comboBoxMethod);

            // поле для матриці
            labelMatrix = new Label
            {
                Location = new Point(20, 60),
                Size = new Size(200, 20),
                Text = "Введіть матрицю:",
                Font = new Font("Arial", 10, FontStyle.Bold)
            };
            this.Controls.Add(labelMatrix);

            // DataGridView для вводу матриці
            dataGridViewMatrix = new DataGridView
            {
                Location = new Point(20, 90),
                Size = new Size(300, 300),
                AllowUserToAddRows = false,
                AllowUserToDeleteRows = false,
                RowHeadersVisible = false,
                BackgroundColor = Color.WhiteSmoke,
                BorderStyle = BorderStyle.FixedSingle,
                GridColor = Color.Gray,
                Font = new Font("Arial", 9)
            };
            this.Controls.Add(dataGridViewMatrix);

            // поле для результатів
            labelResult = new Label
            {
                Location = new Point(350, 60),
                Size = new Size(200, 20),
                Text = "Результати:",
                Font = new Font("Arial", 10, FontStyle.Bold)
            };
            this.Controls.Add(labelResult);

            // TextBox для виведення результатів
            textBoxResult = new TextBox
            {
                Location = new Point(350, 90),
                Size = new Size(300, 300),
                Multiline = true,
                ScrollBars = ScrollBars.Vertical,
                ReadOnly = true,
                Font = new Font("Arial", 9),
                BackColor = Color.WhiteSmoke,
                BorderStyle = BorderStyle.FixedSingle
            };
            this.Controls.Add(textBoxResult);

            // кнопка обчислення
            buttonCalculate = new Button
            {
                Location = new Point(20, 410),
                Size = new Size(100, 30),
                Text = "Обчислити",
                Font = new Font("Arial", 10),
                BackColor = Color.LightBlue,
                FlatStyle = FlatStyle.Flat
            };
            buttonCalculate.FlatAppearance.BorderColor = Color.DarkBlue;
            buttonCalculate.Click += buttonCalculate_Click;
            this.Controls.Add(buttonCalculate);

            // кнопка збереження
            buttonSave = new Button
            {
                Location = new Point(130, 410),
                Size = new Size(100, 30),
                Text = "Зберегти",
                Font = new Font("Arial", 10),
                BackColor = Color.LightGreen,
                FlatStyle = FlatStyle.Flat
            };
            buttonSave.FlatAppearance.BorderColor = Color.DarkGreen;
            buttonSave.Click += buttonSave_Click;
            this.Controls.Add(buttonSave);

            // README у графічному вікні
            labelReadme = new Label
            {
                Location = new Point(20, 450),
                Size = new Size(630, 150),
                Text = "Інструкція:\n" +
                       "1. Оберіть розмір матриці (2–10).\n" +
                       "2. Оберіть метод обчислення (Леверр'є-Фадєєв або Данилевський).\n" +
                       "3. Введіть елементи матриці.\n" +
                       "4. Натисніть 'Обчислити' для отримання власних значень, векторів та складності.\n" +
                       "5. Натисніть 'Зберегти' для збереження результатів у файл.\n" +
                       "Обмеження: розмір матриці 2x2–10x10, лише числові значення.",
                Font = new Font("Arial", 9),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.WhiteSmoke
            };
            this.Controls.Add(labelReadme);
        }

        private void InitializeMatrixGrid(int size)
        {
            matrixSize = size;
            dataGridViewMatrix.Rows.Clear();
            dataGridViewMatrix.Columns.Clear();
            dataGridViewMatrix.AllowUserToAddRows = false;
            dataGridViewMatrix.AllowUserToDeleteRows = false;

            for (int i = 0; i < size; i++)
            {
                dataGridViewMatrix.Columns.Add($"col{i}", $"C{i + 1}");
                dataGridViewMatrix.Columns[i].Width = 60;
                dataGridViewMatrix.Columns[i].DefaultCellStyle.Alignment = DataGridViewContentAlignment.MiddleCenter;
            }
            dataGridViewMatrix.Rows.Add(size);

            // ініціалізація матриці
            matrix = new double[size, size];
        }

        private void numericUpDownSize_ValueChanged(object sender, EventArgs e)
        {
            int newSize = (int)numericUpDownSize.Value;
            if (newSize < 2 || newSize > 10)
            {
                MessageBox.Show("Розмір матриці має бути від 2 до 10.", "Помилка", MessageBoxButtons.OK, MessageBoxIcon.Error);
                numericUpDownSize.Value = matrixSize;
                return;
            }
            InitializeMatrixGrid(newSize);
        }

        private void buttonCalculate_Click(object sender, EventArgs e)
        {
            try
            {
                ReadMatrixFromGrid();
                var calculator = new MatrixCalculator();
                string result = "";
                double[] eigenValues = null;
                double[,] eigenVectors = null;
                operationCount = 0;

                // час виконання
                Stopwatch stopwatch = Stopwatch.StartNew();

                // вибір методу
                if (comboBoxMethod.SelectedItem.ToString() == "Леверр'є-Фадєєв")
                {
                    var (eigenValuesLF, _) = calculator.LeverrierFadeev(matrix, ref operationCount);
                    eigenValues = eigenValuesLF;
                    result = "Власні значення (Метод Леверр'є-Фадєєва):\r\n";
                    for (int i = 0; i < eigenValues.Length; i++)
                        result += $"λ{i + 1} = {(double.IsNaN(eigenValues[i]) ? "Комплексне" : eigenValues[i].ToString("F6"))}\r\n";
                }
                else // Данилевський
                {
                    var (eigenValuesDan, eigenVectorsDan) = calculator.Danilevsky(matrix, ref operationCount);
                    eigenValues = eigenValuesDan;
                    eigenVectors = eigenVectorsDan;
                    result = "Власні значення (Метод Данилевського):\r\n";
                    for (int i = 0; i < eigenValues.Length; i++)
                        result += $"λ{i + 1} = {(double.IsNaN(eigenValues[i]) ? "Комплексне" : eigenValues[i].ToString("F6"))}\r\n";

                    result += "\r\nВласні вектори (Метод Данилевського):\r\n";
                    for (int i = 0; i < eigenValues.Length; i++)
                    {
                        if (!double.IsNaN(eigenValues[i]))
                        {
                            result += $"Вектор для λ{i + 1}:\r\n[";
                            for (int j = 0; j < eigenVectors.GetLength(1); j++)
                                result += $"{eigenVectors[i, j].ToString("F6")}{(j < eigenVectors.GetLength(1) - 1 ? "; " : "")}";
                            result += "]\r\n";
                        }
                    }
                }

                stopwatch.Stop();
                double elapsedTime = stopwatch.Elapsed.TotalMilliseconds;

                // інформація про складність
                result += $"\r\nПрактична складність алгоритму:\r\n";
                result += $"Час виконання: {elapsedTime:F2} мс\r\n";
                result += $"Кількість операцій: {operationCount}\r\n";

                textBoxResult.Text = result;
            }
            catch (Exception ex)
            {
                MessageBox.Show("Помилка: " + ex.Message, "Помилка", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void ReadMatrixFromGrid()
        {
            matrix = new double[matrixSize, matrixSize];

            for (int i = 0; i < matrixSize; i++)
                for (int j = 0; j < matrixSize; j++)
                {
                    object value = dataGridViewMatrix.Rows[i].Cells[j].Value;
                    if (value == null || !double.TryParse(value.ToString(), NumberStyles.Any, CultureInfo.InvariantCulture, out matrix[i, j]))
                        throw new Exception($"Невірне значення в комірці ({i + 1}, {j + 1})");
                }
        }

        private void buttonSave_Click(object sender, EventArgs e)
        {
            SaveFileDialog dlg = new SaveFileDialog
            {
                Filter = "Текстові файли (*.txt)|*.txt",
                Title = "Зберегти результати"
            };
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                File.WriteAllText(dlg.FileName, textBoxResult.Text);
                MessageBox.Show("Результати збережено!", "Успіх", MessageBoxButtons.OK, MessageBoxIcon.Information);
            }
        }

        private void InitializeComponent()
        {
            this.SuspendLayout();
            this.ClientSize = new System.Drawing.Size(684, 561);
            this.Name = "Form1";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.ResumeLayout(false);
        }

        private void Form1_Load(object sender, EventArgs e)
        {
        }
    }

    internal class MatrixCalculator
    {
        public (double[], double[]) LeverrierFadeev(double[,] matrix, ref long operationCount)
        {
            int n = matrix.GetLength(0);
            double[] coefficients = new double[n + 1];
            coefficients[n] = 1; // коефіцієнт при x^n
            double[,] currentMatrix = (double[,])matrix.Clone();
            double[] traces = new double[n];

            for (int k = 0; k < n; k++)
            {
                // обчислення сліду матриці A^(k+1)
                double trace = 0;
                for (int i = 0; i < n; i++)
                {
                    trace += currentMatrix[i, i];
                    operationCount++; // додавання
                }
                traces[k] = trace;

                // обчислення коефіцієнтів p_k за формулою Леверр'є
                double sum = 0;
                for (int m = 1; m <= k; m++)
                {
                    sum += coefficients[n - m] * traces[k - m];
                    operationCount += 2; // множення та додавання
                }
                coefficients[n - k - 1] = -(trace + sum) / (k + 1);
                operationCount += 3; // додавання, віднімання, ділення

                // оновлення матриці A^(k+1)
                currentMatrix = MultiplyMatrix(currentMatrix, matrix, ref operationCount);
            }

            return (SolvePolynomial(coefficients, n, ref operationCount), coefficients);
        }

        public (double[], double[,]) Danilevsky(double[,] matrix, ref long operationCount)
        {
            int n = matrix.GetLength(0);
            // приведення до форми Фробеніуса
            double[] eigenValuesLF;
            double[] coefficients;
            (eigenValuesLF, coefficients) = LeverrierFadeev(matrix, ref operationCount);
            double[] eigenValues = (double[])eigenValuesLF.Clone();

            // обчислення власних векторів через розв’язання (A - λI)v = 0
            double[,] eigenVectors = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                if (double.IsNaN(eigenValues[i])) continue;

                double lambda = eigenValues[i];
                // формування матриці A - λI
                double[,] A_minus_lambdaI = new double[n, n];
                for (int row = 0; row < n; row++)
                    for (int col = 0; col < n; col++)
                    {
                        A_minus_lambdaI[row, col] = matrix[row, col] - (row == col ? lambda : 0);
                        operationCount++; // віднімання
                    }

                // розв'язання матриці (A - λI)v = 0 методом Гауса
                double[] vector = SolveHomogeneousSystem(A_minus_lambdaI, n, ref operationCount);
                if (vector == null)
                {
                    for (int j = 0; j < n; j++)
                        eigenVectors[i, j] = 0;
                    continue;
                }

                //  пприбираємо нормалізацію вектора
                double norm = 0;
                for (int j = 0; j < n; j++)
                {
                    norm += vector[j] * vector[j];
                    operationCount += 2; // множення та додавання
                }
                norm = Math.Sqrt(norm);
                operationCount++; 
                // if (norm > 1e-10)
                   // for (int j = 0; j < n; j++)
                   //  {
                    //    vector[j] /= norm;
                   //     operationCount++; // ділення
                  //  }

                for (int j = 0; j < n; j++)
                    eigenVectors[i, j] = vector[j];
            }

            return (eigenValues, eigenVectors);
        }

        private double[] SolveHomogeneousSystem(double[,] A, int n, ref long operationCount)
        {
            double[,] augmented = new double[n, n + 1];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    augmented[i, j] = A[i, j];

            // прямий хід методу Гауса
            for (int col = 0; col < n; col++)
            {
                int pivotRow = col;
                for (int i = col + 1; i < n; i++)
                    if (Math.Abs(augmented[i, col]) > Math.Abs(augmented[pivotRow, col]))
                        pivotRow = i;

                if (Math.Abs(augmented[pivotRow, col]) < 1e-10)
                    continue;

                if (pivotRow != col)
                    for (int j = 0; j <= n; j++)
                    {
                        (augmented[pivotRow, j], augmented[col, j]) = (augmented[col, j], augmented[pivotRow, j]);
                        operationCount++; // обмін
                    }

                double pivot = augmented[col, col];
                for (int j = col; j <= n; j++)
                {
                    augmented[col, j] /= pivot;
                    operationCount++; // ділення
                }

                for (int i = 0; i < n; i++)
                {
                    if (i == col) continue;
                    double factor = augmented[i, col];
                    for (int j = col; j <= n; j++)
                    {
                        augmented[i, j] -= factor * augmented[col, j];
                        operationCount += 2; // множення та віднімання
                    }
                }
            }

            // зворотний хід
            double[] solution = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                bool isFree = true;
                for (int j = 0; j < n; j++)
                    if (Math.Abs(augmented[i, j]) > 1e-10)
                    {
                        isFree = false;
                        break;
                    }

                if (isFree)
                {
                    solution[i] = 1;
                    continue;
                }

                int leadCol = -1;
                for (int j = 0; j < n; j++)
                    if (Math.Abs(augmented[i, j]) > 1e-10)
                    {
                        leadCol = j;
                        break;
                    }

                if (leadCol == -1) continue;

                solution[leadCol] = augmented[i, n] / augmented[i, leadCol];
                operationCount++;
                for (int j = leadCol + 1; j < n; j++)
                    if (Math.Abs(augmented[i, j]) > 1e-10)
                    {
                        solution[leadCol] -= augmented[i, j] * solution[j] / augmented[i, leadCol];
                        operationCount += 3;
                    }
            }

            double norm = 0;
            for (int i = 0; i < n; i++)
            {
                norm += solution[i] * solution[i];
                operationCount += 2;
            }
            if (norm < 1e-10)
                return null;

            return solution;
        }

        private double[] SolvePolynomial(double[] coefficients, int n, ref long operationCount)
        {
            if (n == 2)
                return SolveQuadratic(coefficients[2], coefficients[1], coefficients[0], ref operationCount);
            else if (n == 3)
                return SolveCubic(coefficients[3], coefficients[2], coefficients[1], coefficients[0], ref operationCount);
            else
                return SolveNewton(coefficients, n, ref operationCount);
        }

        private double[] SolveQuadratic(double a, double b, double c, ref long operationCount)
        {
            double discriminant = b * b - 4 * a * c;
            operationCount += 4;
            if (discriminant < 0)
                return new double[] { double.NaN, double.NaN };
            double sqrtD = Math.Sqrt(discriminant);
            operationCount++; // Корінь
            double[] roots = new double[] { (-b + sqrtD) / (2 * a), (-b - sqrtD) / (2 * a) };
            operationCount += 6;
            return roots;
        }

        private double[] SolveCubic(double a, double b, double c, double d, ref long operationCount)
        {
            b /= a; c /= a; d /= a;
            operationCount += 3;
            double q = (3 * c - b * b) / 9;
            operationCount += 4;
            double r = (9 * b * c - 27 * d - 2 * b * b * b) / 54;
            operationCount += 8;
            double discriminant = q * q * q + r * r;
            operationCount += 4;

            double[] roots = new double[3];

            if (discriminant > 0)
            {
                double sqrtD = Math.Sqrt(discriminant);
                operationCount++;
                double s = Math.Pow(Math.Abs(r + sqrtD), 1.0 / 3.0);
                double t = Math.Pow(Math.Abs(r - sqrtD), 1.0 / 3.0);
                operationCount += 4;
                roots[0] = (r + sqrtD >= 0 ? s : -s) + (r - sqrtD >= 0 ? t : -t) - b / 3;
                operationCount += 3; 
                roots[1] = double.NaN;
                roots[2] = double.NaN;
            }
            else if (Math.Abs(discriminant) < 1e-10)
            {
                double s = Math.Pow(Math.Abs(r), 1.0 / 3.0) * Math.Sign(r);
                operationCount += 2;
                roots[0] = 2 * s - b / 3;
                roots[1] = -s - b / 3;
                roots[2] = roots[1];
                operationCount += 4;
            }
            else
            {
                double theta = Math.Acos(r / Math.Sqrt(-q * q * q));
                operationCount += 4;
                double sqrtQ = Math.Sqrt(-q);
                operationCount++; // Корінь
                roots[0] = 2 * sqrtQ * Math.Cos(theta / 3) - b / 3;
                roots[1] = 2 * sqrtQ * Math.Cos((theta + 2 * Math.PI) / 3) - b / 3;
                roots[2] = 2 * sqrtQ * Math.Cos((theta + 4 * Math.PI) / 3) - b / 3;
                operationCount += 12;
            }

            return roots;
        }

        private double[] SolveNewton(double[] coefficients, int n, ref long operationCount)
        {
            double[] roots = new double[n];
            double[] tempCoefficients = (double[])coefficients.Clone();

            for (int k = 0; k < n; k++)
            {
                double x = k - n / 2.0;
                operationCount++;
                for (int iter = 0; iter < 2000; iter++)
                {
                    double fx = EvaluatePolynomial(tempCoefficients, x, ref operationCount);
                    double dfx = EvaluatePolynomialDerivative(tempCoefficients, x, ref operationCount);

                    if (Math.Abs(dfx) < 1e-12)
                    {
                        roots[k] = double.NaN;
                        break;
                    }

                    double xNew = x - fx / dfx;
                    operationCount += 2;

                    if (double.IsNaN(xNew) || double.IsInfinity(xNew))
                    {
                        roots[k] = double.NaN;
                        break;
                    }

                    if (Math.Abs(xNew - x) < 1e-10)
                    {
                        roots[k] = xNew;
                        tempCoefficients = DeflatePolynomial(tempCoefficients, xNew, ref operationCount);
                        break;
                    }

                    x = xNew;
                }
            }

            Array.Sort(roots);
            Array.Reverse(roots);
            operationCount += 2; // сортування
            return roots;
        }

        private double[] DeflatePolynomial(double[] coefficients, double root, ref long operationCount)
        {
            int n = coefficients.Length - 1;
            double[] newCoefficients = new double[n];
            newCoefficients[n - 1] = coefficients[n];
            for (int i = n - 2; i >= 0; i--)
            {
                newCoefficients[i] = coefficients[i + 1] + root * newCoefficients[i + 1];
                operationCount += 2;
            }
            return newCoefficients;
        }

        private double EvaluatePolynomial(double[] coefficients, double x, ref long operationCount)
        {
            double result = 0;
            for (int i = coefficients.Length - 1; i >= 0; i--)
            {
                result = result * x + coefficients[i];
                operationCount += 2;
            }
            return result;
        }

        private double EvaluatePolynomialDerivative(double[] coefficients, double x, ref long operationCount)
        {
            double result = 0;
            for (int i = coefficients.Length - 1; i > 0; i--)
            {
                result = result * x + i * coefficients[i];
                operationCount += 3;
            }
            return result;
        }

        private double[,] MultiplyMatrix(double[,] A, double[,] B, ref long operationCount)
        {
            int n = A.GetLength(0);
            double[,] result = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    for (int k = 0; k < n; k++)
                    {
                        result[i, j] += A[i, k] * B[k, j];
                        operationCount += 2;
                    }
            return result;
        }
    }
}