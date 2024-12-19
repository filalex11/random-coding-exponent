#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PLOT_DOT_COUNT 10

typedef struct channel_parameters {
	int *ch_in;
	int len_in;
	int *ch_out;
	int len_out;
	double **tr_probs;
} ch_params; 

typedef struct computation_parameters {
	double **phi;
	double *probs_in;
	double **sum;
} comp_params;

/* returns transmit probability of error for pair (y|x) */
double get_tr_prob(double snr, int ch_out, int ch_in) {
	if (ch_in == 1)
		return get_tr_prob(snr, -ch_out, 0);
	double coef = 1.0 / (2.0 * (2.0 + snr));
	if (ch_out > 0)
		return coef * exp(-ch_out / (2.0 * (snr + 1.0)));
	return coef * exp(ch_out / 2.0);
}

void plot_probs() {
	FILE *snr10 = fopen("snr10", "w"),
		*snr25 = fopen("snr25", "w"),
		*snr50 = fopen("snr50", "w"),
		*snr75 = fopen("snr75", "w"),
		*snr100 = fopen("snr100", "w");

	for (double ch_out = -1000; ch_out <= 1000; ch_out += 1) {
		fprintf(snr10, "%e %e\n", ch_out, get_tr_prob(10, ch_out, 0));
		fprintf(snr25, "%e %e\n", ch_out, get_tr_prob(25, ch_out, 0));
		fprintf(snr50, "%e %e\n", ch_out, get_tr_prob(50, ch_out, 0));
		fprintf(snr75, "%e %e\n", ch_out, get_tr_prob(75, ch_out, 0));
		fprintf(snr100, "%e %e\n", ch_out, get_tr_prob(100, ch_out, 0));
	}
}

double find_sum(const ch_params *p, const comp_params *cp) {
	int i, j;
	double res = 0;
	for (i = 0; i < p->len_in; ++i) {
		for (j = 0; j < p->len_out; ++j) {
			res += cp->sum[i][j];
		}
	}
	return res;
}

double get_gal_func(double rho, const ch_params *p, comp_params *cp) {	
	double coef = 1.0 / (1.0 + rho);
	int i, j;
	for (i = 0; i < p->len_in; ++i)
		cp->probs_in[i] = 1.0 / p->len_in;
	for (i = 0; i < p->len_in; ++i) {
		for (j = 0; j < p->len_out; ++j) {
			cp->sum[i][j] = pow(p->tr_probs[i][j], coef) * cp->probs_in[i];
		}
	}
	double res_sum;
	for (i = 0; i < p->len_in; ++i) {
		for (j = 0; j < p->len_out; ++j) {
			res_sum = find_sum(p, cp);
			cp->phi[i][j] = cp->sum[i][j] / res_sum;
			printf("%lf ", cp->phi[i][j]);
		}
		printf("\n");
	}
	for (i = 0; i < p->len_in; ++i) {
		cp->sum[i][0] = 0;
		for (j = 0; j < p->len_out; ++j) {
			cp->sum[i][0] += pow(p->tr_probs[i][j] * 
									pow(cp->phi[i][j], -rho), -1.0 / rho);
		}
	}
	for (i = 0; i < p->len_in; ++i) {
		res_sum = 0;
		for (j = 0; j < p->len_in; ++j)
			res_sum += cp->sum[j][0];
		cp->probs_in[i] = cp->sum[i][0] / res_sum;
		printf("%lf\n", cp->probs_in[i]);
	}
	res_sum = 0;
	for (i = 0; i < p->len_in; ++i) {
		for (j = 0; j < p->len_out; ++j) {
			res_sum += pow(cp->probs_in[i], 1.0 + rho) * p->tr_probs[i][j] 
						* pow(cp->phi[i][j], -rho);
		}
	}
	res_sum = -log(res_sum) * log(2.0);
	return res_sum;
}

comp_params *init_comp_params(const ch_params *p) {
	comp_params *cp = (comp_params*) malloc(sizeof(comp_params));
	cp->phi = (double**) malloc(p->len_in * sizeof(double*));
	cp->sum = (double**) malloc(p->len_in * sizeof(double*));
	cp->probs_in = (double*) malloc(p->len_in * sizeof(double));
	for (int i = 0; i < p->len_out; ++i) {
		cp->phi[i] = (double*) malloc(p->len_out * sizeof(double));
		cp->sum[i] = (double*) malloc(p->len_out * sizeof(double));
	}
	return cp;
}

int main() {
	int i, j, len_in = 2, len_out = 2;
	int ch_in[len_in] = {0, 1};
	int ch_out[len_out];
	double **tr_probs = (double**) malloc(len_in * sizeof(double*));
	double snr = 20;
	for (i = 0; i <= 1; ++i)
		ch_out[i] = i;
	ch_out[0] = -1; ch_out[1] = 1;
	printf("%d %d\n", ch_out[0], ch_out[1]);
	for (i = 0; i < len_in; ++i) {
		tr_probs[i] = (double*) malloc(len_out * sizeof(double));
		for (j = 0; j < len_out; ++j) {
			tr_probs[i][j] = get_tr_prob(snr, ch_out[j], ch_in[i]);
			printf("%lf ", tr_probs[i][j]);
		}
		printf("\n");
	}
	ch_params p;
	p.ch_in = ch_in;
	p.len_in = len_in;
	p.ch_out = ch_out;
	p.len_out = len_out;
	p.tr_probs = tr_probs;
	comp_params *cp = init_comp_params(&p);
	double step = 1.0 / PLOT_DOT_COUNT, rho, value;
	FILE *gal_plot = fopen("gal_func", "w");
	for (i = 0; i <= PLOT_DOT_COUNT; ++i) {
		rho = i * step;
		value = get_gal_func(rho, &p, cp);
		fprintf(gal_plot, "%lf %lf\n", rho, value);
	}
	fclose(gal_plot);
	return 0;
}
