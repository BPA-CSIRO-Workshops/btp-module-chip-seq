class btp {
}

# Helper resource for downloading a remote file
define remote_file($remote_location, $destination, $mode='0644', $owner='root', $group='root') {

  exec { "get_${title}":
    command => "/usr/bin/wget -q ${remote_location} -O ${destination}",
    creates => "${destination}",
    timeout => 0,
  }

  file { "${destination}":
    mode    => $mode,
    owner   => $owner,
    group   => $group,
    require => Exec["get_${title}"],
  }
}

# Helper resource for the analysis tools
define workshop_analysis_tools($source, $provider) {
  remote_file { $name:
    remote_location => "$source/$name",
    destination     => "/usr/local/${name}",
  }

  package { "$source/$name":
    ensure   => present,
    source   => "/usr/local/${name}",
    provider => $provider,
    require  => Remote_file[$name],
  }
}

class btp::analysis_tools {
  $analysis_tools = hiera('btp::analysis_tools', {})

  create_resources(workshop_analysis_tools, $analysis_tools)
}

class btp::system_packages($list=[]) {
  package { $list:
    ensure => present,
  }
}

node image {
  include btp::system_packages
}

node instance {
  include btp::analysis_tools
}
