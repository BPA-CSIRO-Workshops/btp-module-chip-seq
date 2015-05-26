# Helper resource for the directories
define workshop_dir {
  file { "${title}":
    ensure  => directory,
    recurse => true,
    mode    => '0755',
    owner   => $btp::trainee_user,
    group   => $btp::trainee_user,
  }
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

# Helper resource for the workshop files
define workshop_file($location, $link) {
  remote_file { $title:
    remote_location => "${location}/${title}",
    mode            => '0644',
    owner           => $btp::trainee_user,
    group           => $btp::trainee_user,
  }
    
  file { "${link}/${title}":
    ensure  => link,
    target  => "${btp::data_path}/${title}",
    owner   => $btp::trainee_user,
    group   => $btp::trainee_user,
    require => [Workshop_dir[$link], Remote_file[$title]],
  }
}

define workshop_data($source, $targets, $dirs) {
  remote_file { "$data_dir/$name":
    remote_location => "${source}/${name}",
    destination     => "${btp::parent_path}/${btp::data_dir}/${name}",
    owner           => $btp::trainee_user,
    group           => $btp::trainee_user,
    require         => Workshop_dir["${btp::parent_path}/${btp::data_dir}"],
  }

  each($targets) |$target| { 
    remote_file { "$target/$name":
      remote_location => "${source}/${name}",
      destination     => "${btp::parent_path}/${target}",
      owner           => $btp::trainee_user,
      group           => $btp::trainee_user,
      require         => Workshop_subdirs[$dirs],
    }
  }
}

define workshop_subdirs {
  workshop_dir { "${btp::parent_path}/${name}": 
    require => Workshop_dir[$parent_path],
  }
}

define workshop_modules($dirs, $data) {
  $data_path = "${btp::parent_path}/${btp::data_dir}"
  $required_dirs = { 'dirs' => $dirs }

  create_resources(workshop_data, $data, $required_dirs)

  workshop_subdirs { $dirs: 
    before => File["/home/${btp::trainee_user}/${name}"],
  }

  file { [ "/home/${btp::trainee_user}/${name}", 
           "/home/${btp::trainee_user}/Desktop/${name}" ]:
    ensure  => 'link',
    target  => "${btp::parent_path}/${name}",
    owner   => $btp::trainee_user,
    group   => $btp::trainee_user,
    require => Workshop_dir["${btp::parent_path}/${name}"],
  }
}

class btp {
  $parent_path = hiera('btp::parent_path', '/mnt/workshop')
  $data_dir = hiera('btp::data_dir', 'data')
  $trainee_user = hiera('btp::trainee_user', 'trainee')
  $trainee_uid = hiera('btp::trainee_uid', 1001)
  $modules = hiera('btp::modules', {})

  $data_path = "${parent_path}/${data_dir}"

  # Defaults
  File {
    owner => $trainee_user,
    group => $trainee_user,
  }

  # Trainee group
  group { $trainee_user:
    ensure => present,
  }

  # Trainee user's home directory
  workshop_dir { "/home/${trainee_user}": }

  # Trainee user's desktop directory
  workshop_dir { "/home/${trainee_user}/Desktop":
    require => Workshop_dir["/home/${trainee_user}"],
  }

  # Trainee user
  user { $trainee_user:
    ensure  => present,
    uid     => $trainee_uid,
    gid     => $trainee_user,
    shell   => '/bin/bash',
    home    => "/home/${trainee_user}",
    require => Group[$trainee_user],
  }

  # Parent path
  workshop_dir { $parent_path:
    require => User[$trainee_user],
  }

  # Data directory
  workshop_dir { $data_path:
    require => Workshop_dir[$parent_path],
  }

  create_resources(workshop_modules, $modules)
}

node default {
  include btp
}
