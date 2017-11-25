OPTIONS_AUTOMAKE=gnu
ACLOCAL_AMFLAGS=-I ./m4

bin_PROGRAMS = anova correct fit glm ks nna pca svd ttest

lib_LIBRARIES = libxa_stats.a

libxa_stats_a_SOURCES = files.c function.c matrix.c methods.c metrics.c solution.c stats.c utils.c var_metrics.c ttest_impl.c
libxa_stats_a_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
libxa_stats_a_includedir=${includedir}/xa_stats
libxa_stats_a_include_HEADERS = metrics.h solution.h utils.h function.h matrix.h stats.h var_metrics.h files.h methods.h ttest_impl.h xa_stats_version.h
libxa_stats_a_DEPENDENCIES =

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = xa_stats.pc
pkginclude_HEADERS = ${libxa_stats_a_include_HEADERS}
pkginclude_LIBS = libxa_stats

anova_SOURCES = anova.c
anova_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
anova_LDFLAGS =
anova_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
anova_DEPENDENCIES = libxa_stats.a anova.h

correct_SOURCES = correct.c
correct_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
correct_LDFLAGS =
correct_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
correct_DEPENDENCIES = libxa_stats.a correct.h ${libxa_stats_a_include_HEADERS}

fit_SOURCES = fit.c
fit_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
fit_LDFLAGS =
fit_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
fit_DEPENDENCIES = libxa_stats.a fit.h ${libxa_stats_a_include_HEADERS}

glm_SOURCES = glm.c
glm_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
glm_LDFLAGS =
glm_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
glm_DEPENDENCIES = libxa_stats.a glm.h ${libxa_stats_a_include_HEADERS}

ks_SOURCES = ks.c
ks_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
ks_LDFLAGS =
ks_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
ks_DEPENDENCIES = libxa_stats.a ks.h ${libxa_stats_a_include_HEADERS}

nna_SOURCES = nna.c
nna_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
nna_LDFLAGS =
nna_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
nna_DEPENDENCIES = libxa_stats.a nna.h ${libxa_stats_a_include_HEADERS}

pca_SOURCES = pca.c
pca_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
pca_LDFLAGS =
pca_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
pca_DEPENDENCIES = libxa_stats.a pca.h ${libxa_stats_a_include_HEADERS}

svd_SOURCES = svd.c
svd_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
svd_LDFLAGS =
svd_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
svd_DEPENDENCIES = libxa_stats.a svd.h ${libxa_stats_a_include_HEADERS}

ttest_SOURCES = ttest.c
ttest_CFLAGS = -I. ${OPTIMISATION_FLAGS} ${CFLAGS}
ttest_LDFLAGS =
ttest_LDADD = -L. ${LDFLAGS} -lxa_stats -lm
ttest_DEPENDENCIES = libxa_stats.a ttest.h ${libxa_stats_a_include_HEADERS}
