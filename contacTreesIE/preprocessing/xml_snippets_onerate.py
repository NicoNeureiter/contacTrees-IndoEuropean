

# Alignments
ALIGNMENT = '''\
    <sequence id=\"language_data:{tax_name}\" taxon=\"{tax_name}\" value=\"{data}\"/>
'''
FILTERED_ALIGNMENT = '''\
  <data id="feature_data:{concept}" spec="FilteredAlignment" data="@data" filter="{start}-{end}" ascertained="true" excludeto="{excludeto}" userDataType="@dtype"/>
'''

# Likelihood

BASICTREES_LIKELIHOOD = '''\
          <distribution id="treeLikelihood.$(concept)" spec="TreeLikelihood" useAmbiguities="true">
            <tree idref="acg" />
            <branchRateModel idref="clock" />
            <siteModel idref="SiteModel.vocabulary" />
            <data idref="feature_data:$(concept)" />
          </distribution>
'''

CONTACTREES_LIKELIHOOD = '''\
          <distribution id="treeLikelihood.$(concept)" spec="TreeLikelihood" useAmbiguities="true">
            <tree id="marginalTree.$(concept)" spec="MarginalTree" network="@acg" block="@$(concept)" branchRateModel="@clock" />
            <siteModel idref="SiteModel.vocabulary" />
            <data idref="feature_data:$(concept)" />
          </distribution>
'''


# Operators

TOPOLOGY_OPERATORS = '''\
    <operator id="CFWilsonBalding" spec="CFWilsonBalding" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" includeRoot="false" alpha="0.1" networkPrior="@ACGBirthDeathModel" weight="1.0"/>
    <operator id="CFNarrowExchange" spec="CFSubtreeExchange" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" isNarrow="true" networkPrior="@ACGBirthDeathModel" weight="10.0"/>
    <operator id="CFWideExchange" spec="CFSubtreeExchange" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" isNarrow="false" networkPrior="@ACGBirthDeathModel" weight="1.0"/>
'''

NODE_HEIGHTS_OPERATORS = '''\
    <operator id="ACGscaler" spec="ACGScaler" acg="@acg" scaleFactor="0.9" weight="3.0"/>
    <operator id="ACGscaler.rootOnly" spec="ACGScaler" acg="@acg" scaleFactor="0.85" weight="1.0" rootOnly="true"/>
    <operator id="CFUniform" spec="CFUniform" acg="@acg" conversionRate="@conversionRate" pMove="@pMove" blockSet="@allBlocks" scaleFactor="0.9" networkPrior="@ACGBirthDeathModel" weight="28.0"/>
'''

CONTACT_OPERATORS = '''\
    <!--operator id="AddRemoveConversion.t" spec="AddRemoveConversion" weight="50.0" acg="@acg" pMove="@pMove" conversionRate="@conversionRate" blockSet="@allBlocks" networkPrior="@ACGBirthDeathModel"/-->
    <operator id="AddRemoveConversion.t" spec="AddRemoveConversionGibbs" weight="50.0" acg="@acg" pMove="@pMove" conversionRate="@conversionRate" blockSet="@allBlocks" networkPrior="@ACGBirthDeathModel">
        <plate var="concept" range="&concepts;">
          <treeLikelihood idref="treeLikelihood.$(concept)"/>
        </plate>
    </operator>
    <operator id="GibbsSampleMovesPerConversion.t" spec="GibbsSampleMovesPerConversion" weight="10.0" acg="@acg" pMove="@pMove" blockSet="@allBlocks">
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
    <operator id="ConvertedEdgeHop.source.narrow" spec="ConvertedEdgeHopGibbs" acg="@acg" sourceOnly="true" nClosestRelatives="3" blockSet="@allBlocks" pMove="@pMove" conversionRate="@conversionRate" networkPrior="@ACGBirthDeathModel" weight="5.0">
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

# Loggers

WORD_TREE_LOGGERS = '''
    <plate var="concept" range="&concepts;">
      <logger spec="Logger" fileName="wordtrees/$(filebase).$(concept).trees" logEvery="125000" mode="tree">
        <log spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@marginalTree.$(concept)"/>
      </logger>
    </plate>
'''
