$:.unshift File.join(File.dirname(__FILE__),'..','lib')

require 'rubygems'
require 'test/unit'
require 'go'

class GoTest < Test::Unit::TestCase
  def setup
    @go = Bio::Go.new
  end
  
  def test_cellular_component_offspring
    # test no offspring
    # GO:0031676 does not (as of test writing) have any offspring
    assert_equal [], @go.cellular_component_offspring('GO:0031676')
    
    # test multiple offspring
    old_offspring = ["GO:0030077", "GO:0030078", "GO:0030079", "GO:0030080", "GO:0030081", "GO:0030082", "GO:0031633", "GO:0031676", "GO:0042717", "GO:0048493", "GO:0048494"]
    offspring = @go.cellular_component_offspring('GO:0042716')
    assert offspring.length >= old_offspring.length
    assert offspring.include?(old_offspring[0])
    assert offspring.include?(old_offspring[1])
    assert offspring.include?(old_offspring[old_offspring.length-2])
    
    # test not in CC
    assert_raise RException do
      @go.cellular_component_offspring('GO:0042716not')
    end
  end
  
  # operations below are expensive because they require
  # multiple downloads from internet, so are commented
  # out by default
  #  def test_cc_pdb_to_go_some
  #    assert_equal ['GO:0005743'],
  #      @go.cc_pdb_to_go('2a06')
  #  end
  
  def test_go_term
    # test MF
    assert_equal "G-protein coupled receptor activity", @go.term('GO:0004930')
    
    # test CC
    assert_equal 'endoplasmic reticulum', @go.term('GO:0005783')
  end
  
  def test_synonym
    # test real synonym
    assert_equal 'GO:0050333', @go.primary_go_id('GO:0048253')
    
    # test primary id
    assert_equal 'GO:0050333', @go.primary_go_id('GO:0050333')
    
    # test bad id
    assert_raise RException do
      @go.primary_go_id('GO:AWAY')
    end
  end
  
  def test_subsume
    # test normal truth
    assert @go.subsume?('GO:0003824', 'GO:0050333')
    
    # test subsumee is synonym
    assert @go.subsume?('GO:0003824', 'GO:0048253')
    
    # test equal terms
    assert @go.subsume?('GO:0009536','GO:0009536')
    
    # test falsity - plastid part does not subsume plastid
    assert_equal false, @go.subsume?('GO:0044435','GO:0009536')
  end
  
  def test_subsume_tester
    tester = @go.subsume_tester('GO:0003824')
    assert_kind_of Bio::Go::SubsumeTester, tester
    assert tester.subsume?('GO:0050333')
    
    tester = @go.subsume_tester('GO:0044435')
    assert_equal false, tester.subsume?('GO:0009536')
    
    #equality
    tester = @go.subsume_tester('GO:0044435')
    assert tester.subsume?('GO:0044435')
  end
  
  def test_subsume_tester_no_check_synonym
    # test normal truth
    tester = @go.subsume_tester('GO:0003824')
    assert_equal true, tester.subsume?('GO:0050333', false)
    assert_equal true, tester.subsume?('GO:0050333', true)
    
    # test subsumee is synonym
    tester = @go.subsume_tester('GO:0003824')
    assert_equal false, tester.subsume?('GO:0048253', false)
    assert_equal true, tester.subsume?('GO:0048253', true)
    
    # test equal terms
    tester = @go.subsume_tester('GO:0050333')
    assert_equal true, tester.subsume?('GO:0050333', false)
    assert_equal true, tester.subsume?('GO:0050333', true)
    
    # test equal that is synonym
    tester = @go.subsume_tester('GO:0050333')
    assert_equal false, tester.subsume?('GO:0048253', false)
    assert_equal true, tester.subsume?('GO:0048253', true)
  end
  
  def test_cordial_cc
    nucleus_cordial = @go.cordial_cc('GO:0005634')
    
    assert nucleus_cordial.include?('GO:0005575'), "include cellular component itself"
    assert_equal false,  nucleus_cordial.include?('GO:0005634'), 'don\'t include nucleus itself'
    assert nucleus_cordial.include?('GO:0005623'), 'include a direct ancestor of nucleus'
    assert nucleus_cordial.include?('GO:0005623'), 'include membrane, a \'indirect\' direct ancestor of nucleus (it is an ancestor of nuclear membrane)'
  end
end
