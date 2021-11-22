function pval = r2_barplot(r2_output)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

figure

sumstats_expected_se = r2_output.sumstats_expected_se;
mtag_expected_se = r2_output.mtag_expected_se;
mtag_observed_se = r2_output.mtag_observed_se;
expected_se = r2_output.pdr_expected_se;
observed_se = r2_output.pdr_observed_se;
sumstats_expected = r2_output.sumstats_expected;
mtag_expected = r2_output.mtag_expected;
mtag_observed = r2_output.mtag_observed;
expected = r2_output.pdr_expected;
observed = r2_output.pdr_observed;
pdr_observed_jk = r2_output.pdr_observed_jk;
mtag_observed_jk = r2_output.mtag_observed_jk;

if isempty(r2_output.pdr_observed)
    mtag_observed = 0;
    mtag_observed_se = 0;
    observed = 0;
    observed_se = 0;
end

barwitherr([sumstats_expected_se, 0;
    mtag_expected_se mtag_observed_se;
    expected_se observed_se],...
    [sumstats_expected, 0;
    mtag_expected mtag_observed;
    expected observed])
ylim([0 1])
ylabel('Variance explained (r^2)')
legend('Expected r^2','Observed r^2')
set(gca,'XTickLabel',{'Sumstats','MTAG','PDR'})
pval = 2*normcdf(-mean(pdr_observed_jk-mtag_observed_jk)/std(pdr_observed_jk-mtag_observed_jk)/10);


end

