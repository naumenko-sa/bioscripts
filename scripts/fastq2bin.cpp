#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

//fastq2binary - by T.Pirogova

void compress(FILE *inf, FILE *outf)
{
	unsigned char buf[1000], str[255];
	unsigned char ch = 0, s;
	std::string errStr; 
	fwrite(&ch, 1, 1, outf);
	bool isErr, first = true;
	while(!feof(inf))
	{
		ch = 0;
		isErr = false;
		fgets((char*)buf, 1000, inf);
		if (first) errStr.clear();
		errStr += (char*)buf;
                if (feof(inf)) return;
		fgets((char*)buf, 1000, inf);
		errStr += (char*)buf;
		unsigned short length =  strlen((char*)buf)-1;
		buf[length] = 0;
		unsigned char *tmp = str;
		for (unsigned short i = 0; i < length; i++)
		{
			switch (buf[i])
			{
				case 'A': 
				case 'a': s = 0;
					break;
				case 'C':
				case 'c': s = 1;
					break;
				case 'G':
				case 'g': s = 2;
					break;
				case 'T':
				case 't': s = 3;
					break;
				default: isErr = true;
			}
			if (isErr) break;
			ch |= (s << ((i%4)<<1) ); // добавляем следующую букву, смещая ее на 0,2,4,6 бит в зависимости от i
			if ( (!((i+1)%4)) || ((i+1) == length) ) // если 4-й элемент или последний
			{
				*tmp++ = ch; // в противном случае добавляем в массив готовый байт
				ch = 0;
			}
		}
		fgets((char*)buf, 1000, inf);
		errStr += (char*)buf;
		fgets((char*)buf, 1000, inf);
		errStr += (char*)buf;
		if (isErr) 
		{
			if (first)
				for (int i = 0; i<4; i++)
				{
					fgets((char*)buf, 1000, inf);
					errStr += (char*)buf;
				}
			printf("%s\n", errStr.c_str());
			continue;
		}
		first = !first;
		fwrite(&length, sizeof(short), 1, outf);
		fwrite(str, (length+3)/4, 1, outf);
	}
}

void decompress(FILE *inf, FILE *outf)
{
	unsigned short length, count, i;
	unsigned char buf[1000], s, str[1000];
	bool isErr = false;
	while (!feof(inf))
	{
		fread(&length, sizeof(short), 1, inf);
		fread(buf, (length+3)/4, 1, inf);
		for (i = 0; i<length;)
		{
			s = buf[i/4];
			for (count = 0; count<4; count++, i++)
			{
				if (i>=length) break;			
				switch (s&3)
				{
				case 0:	str[i] = 'A';
					break;
				case 1: str[i] = 'C';
					break;
				case 2: str[i] = 'G';
					break;
				case 3: str[i] = 'T';
					break; 
				default: isErr = true;
				}
				if (isErr)
				{
					printf("Unresolved simbol in decomressing file\n");
					return;
				}
				s = s>>2;
			}
		}
		str[i] = '\n';
		str[++i] = 0;
		fputs((char*)str, outf);
	}
}

int main (int argc, char *argv[])
{
	std::string fname;
	FILE *inf, *outf;
	bool isOk = true;
	if (argc < 2)
	{
		printf("Input file is not defined.\n""%s input_file\n", argv[0]);
		return 0;
	}
	else fname = argv[1];
	if ((inf = fopen(fname.c_str(),"rb")) == NULL)
	{
            printf("The file %s was not opened\n", fname.c_str());
		return 0;
	}
	char buf;
	fread(&buf,1,1,inf);
	if (buf == '\0') 
	{
		fname += ".decompressed";
		if ((outf = fopen(fname.c_str(), "wt")) == NULL)
		{
                    printf("The file %s was not opened\n", fname.c_str());
			isOk = false;
		}
		else 
			decompress(inf, outf);
	}
	else
	{
		fclose(inf);
		if ((inf = fopen(fname.c_str(),"rt")) == NULL)
		{
                    printf("The file %s was not opened\n", fname.c_str());
			return 0;
		}
		if ((outf = fopen((fname+".compressed").c_str(), "wb")) == NULL)
		{
                    printf("The file %s was not opened\n", fname.c_str());
			isOk = false;
		}
		else
			compress(inf, outf);
	}
	fclose(inf);
	if (isOk)
                fclose(outf);
	return 0;
}

