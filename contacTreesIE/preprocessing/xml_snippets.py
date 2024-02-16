from enum import Enum

# Alignments
ALIGNMENT = '''\
  <sequence id=\"language_data:{tax_name}\" taxon=\"{tax_name}\" value=\"{data}\"/>
'''
FILTERED_ALIGNMENT = '''\
<data id="feature_data:{concept}" spec="FilteredAlignment" data="@data" filter="{start}-{end}" ascertained="true" excludeto="{excludeto}" userDataType="@dtype"/>
'''


class Samplers(Enum):
    MCMC = 1
    MC3 = 2
    NS = 3

SAMPLER_TAG = {
    Samplers.MCMC: '<run id="mcmc" spec="MCMC" chainLength="{chain_length}">',
    Samplers.MC3: '<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" chainLength="{chain_length}" chains="4" deltaTemperature="0.04" target="0.15" optimise="true" resampleEvery="2000" heatLikelihoodOnly="true" logHeatedChains="true" neighbourSwapping="true">',
    Samplers.NS: '<run id="mcmc" spec="beast.gss.MultiThreadedNS" chainLength="{chain_length}" threads="2" particleCount="1" subChainLength="50000">',
    # Samplers.NS: '<run id="mcmc" spec="beast.gss.NS" chainLength="{chain_length}" particleCount="1" subChainLength="20000">',
}
SAMPLER_THREADS = {
    Samplers.MCMC: 1,
    Samplers.MC3: 4,
    Samplers.NS: 2,
}


# Likelihood

BASICTREES_LIKELIHOOD = '''\
      <distribution id="treeLikelihood.{concept}" spec="TreeLikelihood" useAmbiguities="true" tree="@acg" 
          branchRateModel="@clock.{site_cat}" siteModel="@SiteModel.vocabulary" data="@feature_data:{concept}"/>
'''

CONTACTREES_LIKELIHOOD = '''\
      <distribution id="treeLikelihood.{concept}" spec="TreeLikelihood" useAmbiguities="true" siteModel="@SiteModel.vocabulary" data="@feature_data:{concept}">
        <tree id="marginalTree.{concept}" spec="MarginalTree" network="@acg" block="@{concept}" branchRateModel="@clock.{site_cat}"/>
      </distribution>
'''
TMP = '''\
        <tree id="marginalTree.{concept}" spec="MarginalTree" network="@acg" block="@{concept}" branchRateModel="@clock.{site_cat}">
            <frozenTaxa>Latin_preserved</frozenTaxa>
        </tree>
'''

FROZEN_TAXA = '  <frozenTaxa idref="Latin_preserved"/>'

# Clock model
CLOCK_RATE_PRIOR_FLAT = '''\
      <prior id="prior:clockRate.slow" name="distribution" x="@clockRate.slow">
        <distr id="prior:clockRate.slow:distr" spec="beast.math.distributions.Uniform" lower="0.0" upper="1.0" />
      </prior>
      <prior id="prior:clockRate.medium" name="distribution" x="@clockRate.medium">
        <distr id="prior:clockRate.medium:distr" spec="beast.math.distributions.Uniform" lower="0.0" upper="1.0" />
      </prior>
      <prior id="prior:clockRate.fast" name="distribution" x="@clockRate.fast">
        <distr id="prior:clockRate.fast:distr" spec="beast.math.distributions.Uniform" lower="0.0" upper="1.0" />
      </prior>\
'''
CLOCK_RATE_PRIOR_INFOMRATIVE = '''\
      <prior id="prior:clockRate.slow" name="distribution">
        <x id="lossRate.slow" spec="feast.expressions.ExpCalculator">
           <arg idref="clockRate.slow"/>
           <arg idref="freqParameter"/>
           clockRate.slow / (2 * freqParameter[1])
        </x>
        <distr id="prior:clockRate.slow:distr" spec="LogNormalDistributionModel" offset="0.0" M="-8.5" S="2.3" />
      </prior>
      <prior id="prior:clockRate.medium" name="distribution">
        <x id="lossRate.medium" spec="feast.expressions.ExpCalculator">
           <arg idref="clockRate.medium"/>
           <arg idref="freqParameter"/>
           clockRate.medium / (2 * freqParameter[1])
        </x>
        <distr id="prior:clockRate.medium:distr" spec="LogNormalDistributionModel" offset="0.0" M="-7.5" S="2.3" />
      </prior>
      <prior id="prior:clockRate.fast" name="distribution">
        <x id="lossRate.fast" spec="feast.expressions.ExpCalculator">
           <arg idref="clockRate.fast"/>
           <arg idref="freqParameter"/>
           clockRate.fast / (2 * freqParameter[1])
        </x>
        <distr id="prior:clockRate.fast:distr" spec="LogNormalDistributionModel" offset="0.0" M="-6.5" S="2.3" />
      </prior>\
'''

# Substitution model

CTMC_MODEL = '''\
<substModel id="substModel" spec="GeneralSubstitutionModel">
  <rates spec="RealParameter" dimension="2" estimate="false" id="rates.binaryCTMC">1.0 1.0</rates>
  <frequencies id="frequenciesCTMC" spec="Frequencies" frequencies="@freqParameter"/>
</substModel>\
'''
CTMC_PRIORS = ''
CTMC_DATA_TYPE = '<userDataType id="dtype" spec="beast.evolution.datatype.Binary" />'
COVARION_MODEL = '''\
<substModel id="substModel" spec="BinaryCovarion" alpha="@covarionAlpha" switchRate="@covarionSwitchRate" vfrequencies="@freqParameter">
  <parameter id="hiddenfrequencies" dimension="2" lower="0.0" name="hfrequencies" upper="1.0">0.5 0.5</parameter>
  <frequencies id="dummyfrequencies" spec="Frequencies" frequencies="0.5 0.5" />
</substModel>\
'''
COVARION_PRIORS = '''\
      <prior id="prior:covarionAlpha" name="distribution" x="@covarionAlpha">
        <distr id="prior:covarionAlpha:distr" spec="beast.math.distributions.Uniform" lower="0" upper="1" />
      </prior>

      <prior id="prior:covarionSwitchRate" name="distribution" x="@covarionSwitchRate">
        <distr id="prior:covarionSwitchRate:distr" spec="Exponential" mean="0.2" />
      </prior>
'''
COVARION_DATA_TYPE = '<userDataType id="dtype" spec="beast.evolution.datatype.TwoStateCovarion" />'


# Operators

TOPOLOGY_OPERATORS = '''\
    <operator id="CFWilsonBalding" spec="CFWilsonBalding" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" includeRoot="false" alpha="0.1" networkPrior="@ACGBirthDeathModel" weight="1.0"/>
    <operator id="CFNarrowExchange" spec="CFSubtreeExchange" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" isNarrow="true" networkPrior="@ACGBirthDeathModel" weight="10.0"/>
    <operator id="CFWideExchange" spec="CFSubtreeExchange" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" isNarrow="false" networkPrior="@ACGBirthDeathModel" weight="1.0"/>
'''

NODE_HEIGHTS_OPERATORS = '''\
    <operator id="ACGscaler" spec="ACGScaler" acg="@acg" scaleFactor="0.9" weight="10.0">
        <parameterInverse idref="clockRate.slow"/>
        <parameterInverse idref="clockRate.medium"/>
        <parameterInverse idref="clockRate.fast"/>
    </operator>
    <operator id="ACGscaler.rootOnly" spec="ACGScaler" acg="@acg" scaleFactor="0.75" weight="1.0" rootOnly="true"/>
    <operator id="CFUniform" spec="CFUniform" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" scaleFactor="0.9" networkPrior="@ACGBirthDeathModel" weight="28.0"/>
'''

CONTACT_OPERATORS = '''\
    <!--operator id="AddRemoveConversion.t" spec="AddRemoveConversion" weight="50.0" acg="@acg" pMove="@pMove" conversionRate="@conversionRate" blockSet="@allBlocks" networkPrior="@ACGBirthDeathModel"/-->
    <operator id="AddRemoveConversion.t" spec="AddRemoveConversionGibbs" weight="50.0" acg="@acg" pMove="@pMove" conversionRate="@conversionRate" blockSet="@allBlocks" networkPrior="@ACGBirthDeathModel">
        <plate var="concept" range="&concepts;">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="GibbsSampleMovesPerConversion.t" spec="GibbsSampleMovesPerConversion" weight="10.0" acg="@acg" pMove="@pMove" blockSet="@allBlocks" mcmcmc="true">
        <plate var="concept" range="&concepts;">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>

    <operator id="ConvertedEdgeSlide.t" spec="ConvertedEdgeSlide" acg="@acg" apertureSize="0.2" weight="15.0"/>
    <operator id="ConvertedEdgeFlip.t" spec="ConvertedEdgeFlip" acg="@acg" weight="1.0"/>
    <operator id="ConversionSplit.t" spec="ConversionSplit" acg="@acg" weight="1.0"
            blockSet="@allBlocks" conversionRate="@conversionRate" networkPrior="@ACGBirthDeathModel" flip="false" pMove="@pMove"/>
    <operator id="ConvertedEdgeHop.source" spec="ConvertedEdgeHopGibbs" acg="@acg" sourceOnly="true" blockSet="@allBlocks" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGBirthDeathModel" weight="2.0">
        <plate var="concept" range="&concepts;">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="ConvertedEdgeHop.source.narrow" spec="ConvertedEdgeHopGibbs" acg="@acg" sourceOnly="true" nClosestRelatives="3" blockSet="@allBlocks" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGBirthDeathModel" weight="6.0">
        <plate var="concept" range="&concepts;">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="ConvertedEdgeHop.narrow" spec="ConvertedEdgeHopGibbs" acg="@acg" blockSet="@allBlocks" nClosestRelatives="3" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGBirthDeathModel" weight="2.0">
        <plate var="concept" range="&concepts;">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
'''

CLOCK_STDEV_OPERATOR = '  <operator id="ClockStdevScaler" spec="ScaleOperator" parameter="@clockStdev" scaleFactor="0.7" weight="4.0"/>'

# Loggers

WORD_TREE_LOGGERS = '''
    <plate var="concept" range="&concepts;">
      <logger spec="Logger" fileName="wordtrees/$(filebase).$(concept).trees" logEvery="125000" mode="tree">
        <log spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@marginalTree.$(concept)"/>
      </logger>
    </plate>
'''
