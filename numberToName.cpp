const string X[] = {"", "One ", "Two ", "Three ", "Four ", "Five ",
                    "Six ", "Seven ", "Eight ", "Nine ", "Ten ", "Eleven ",
                    "Twelve ", "Thirteen ", "Fourteen ", "Fifteen ",
                    "Sixteen ", "Seventeen ", "Eighteen ", "Nineteen "};

const string Y[] = {"", "", "Twenty ", "Thirty ", "Forty ", "Fifty ",
                    "Sixty ", "Seventy ", "Eighty ", "Ninety "};

string convertToDigit(int n, string suff)
{
    if (n == 0)
    {
        return "";
    }

    if (n > 19)
    {
        return Y[n / 10] + X[n % 10] + suff;
    }
    else
    {
        return X[n] + suff;
    }
}

string nameTheNumber(int n)
{
    if (n == 0)
    {
        return "Zero";
    }

    if (n == 1000000000)
    {
        return "One Billion";
    }

    string res;

    res = convertToDigit((n % 100), "");

    res = convertToDigit(((n / 100) % 10), "Hundred ") + res;

    res = convertToDigit(((n / 1000) % 100), "Thousand ") + res;

    res = convertToDigit(((n / 100000) % 10), "Hundred ") + res;

    res = convertToDigit(((n / 1000000) % 100), "Million ") + res;

    res = convertToDigit(((n / 100000000) % 10), "Hundred ") + res;

    res = convertToDigit(((n / 1000000000) % 100), "Billion ") + res;

    return res;
}